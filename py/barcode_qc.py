from argparse import ArgumentParser
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import seaborn as sbn
import numpy as np
import gzip

def get_args():
    parser = ArgumentParser()
    parser.add_argument('annot_fq', help='The annotated fastq file', nargs='+')
    parser.add_argument('--output_pdf', help='The output qc pdf', default=None)
    parser.add_argument('--plate_layout', help='The plate layout csv file', default='/home/chartl/repos/pipelines/config/plate_layout.csv')

    return parser.parse_args()


def main(args):
    libname = args.annot_fq[1].split('.fq')[0].split('/')[-1].split('__')[0]
    if args.output_pdf is None:
        outpdf = libname + '.barcode_qc.pdf'
    else:
        outpdf = args.output_pdf

    def title(ttl):
        plt.title(ttl + '\nLibrary: %s' % libname)

    counts = dict()
    for fastq in args.annot_fq:
        hdl = gzip.open(fastq, 'rt')
        for i, line in enumerate(hdl):
            if (i % 4) != 0:
                continue
            vals = line.strip().split('|')[1].split(':')
            if vals[0] == '*':
                key = ('*', '*')
            elif '*' in vals[1] or '*' in vals[2]:
                key = (vals[3], '*')
            else:
                key = (vals[3], vals[1] + '.' + vals[2])
            counts[key] = 1 + counts.get(key, 0)
        hdl.close()

    records = pd.DataFrame([{
       'sample': k[0],
       'cell': k[1],
       'sample_cell': k[0] + ':' + k[1],
       'reads': v} for k, v in counts.items()])

    # plate heatmaps
    plate_table = pd.read_csv(args.plate_layout)
    barcode2info = {(row.barcode, row.split): (row.plate, row.plate_row, row.plate_col) for ix_, row in plate_table.iterrows()}

    def get_info(bc_id, split, what='plate'):
        if bc_id == '*':
            return 999
        if (bc_id, split) in barcode2info:
            plate, row, col = barcode2info[(bc_id, split)]
        else:
            return -1
        if what == 'plate':
            return int(plate[1:])
        elif what == 'row':
            return 1 + row
        else:
            return 1 + col
    
    records.loc[:, 'well1_plate'] = records.cell.map(lambda x: get_info(x.split('.')[0], 1, what='plate'))
    records.loc[:, 'well2_plate'] = records.cell.map(lambda x: get_info(x.split('.')[1] if '.' in x else '*', 2, what='plate'))
    records.loc[:, 'well1_row'] = records.cell.map(lambda x: get_info(x.split('.')[0], 1, what='row'))
    records.loc[:, 'well2_row'] = records.cell.map(lambda x: get_info(x.split('.')[1] if '.' in x else '*', 2, what='row'))
    records.loc[:, 'well1_col'] = records.cell.map(lambda x: get_info(x.split('.')[0], 1, what='col'))
    records.loc[:, 'well2_col'] = records.cell.map(lambda x: get_info(x.split('.')[1] if '.' in x else '*', 2, what='col'))
    records.loc[:, 'assignment_valid'] = ((records.well1_plate != -1) & (records.well2_plate != -1)).map(lambda x: 'valid' if x else 'invalid')

    r_unassign = records[records.sample_cell.map(lambda x: '*' in x)].reads.sum() / records.reads.sum()
    str_A = 'assigned (%.1f%%)' % (100 * (1-r_unassign)) 
    str_U = 'unassigned (%.1f%%)' % (100 * r_unassign)
    records.loc[:, 'assignment'] = records.sample_cell.map(lambda x: str_U if '*' in x else str_A)
    acounts = records.groupby('assignment')['reads'].sum().reset_index()
    pdf = PdfPages(outpdf)
    plt.bar(acounts.assignment, acounts.reads)
    plt.ticklabel_format(axis='y', style='plain')

    plt.tight_layout()
    pdf.savefig()
    plt.figure()

    n_valid = records[(records.well1_plate != 999) & (records.assignment_valid == 'valid')].shape[0]
    n_invalid = records[(records.well1_plate != 999) & (records.assignment_valid == 'invalid')].shape[0]

    valid_counts = records[records.sample_cell.map(lambda x: '*' not in x)].assignment_valid.value_counts()
    plt.bar(['valid barcodes', 'invalid barcodes'], [n_valid, n_invalid])
    plt.tight_layout()
    pdf.savefig()
    plt.figure()

    library_total_reads = acounts[acounts.assignment == str_A].reads

    scounts = records.groupby(['sample', 'assignment'])['reads'].sum().reset_index()
    sbn.barplot(x='sample', y='reads', hue='assignment', data=scounts[scounts['sample'] != '*'])
    plt.xticks(rotation=90)
    plt.ticklabel_format(axis='y', style='plain')
    plt.tick_params(labelsize=8)
    ax = plt.gca()
    ax.get_legend().remove()

    plt.tight_layout()
    pdf.savefig()
    plt.figure()

    kmin, kmax, kstep = 250, 2000, 50
    dat = records[records.assignment == str_A]
    dat_cdf = pd.DataFrame(dict(n_reads=np.arange(kmin, kmax, kstep),
                    n_cells=[(dat.reads >= x).sum() 
                    for x in np.arange(kmin, kmax, kstep)])) 
    plt.plot(dat_cdf.n_reads, dat_cdf.n_cells)
    plt.plot([500, 500], [dat_cdf.n_cells.min() * 0.8, dat_cdf.n_cells.max() * 1.05], color='black')
    plt.ylabel('# Cells with >X reads')
    plt.xlabel('# Reads assigned to cell')

    plt.tight_layout()
    pdf.savefig()
    plt.figure()

    dat_cdf_persample = {'n_reads': np.arange(kmin, kmax, kstep)}
    for sn in dat['sample'].unique():
        dat_cdf_persample[sn] = [(dat[dat['sample'] == sn].reads > x).sum()
                         for x in np.arange(kmin, kmax, kstep)]
        plt.plot(dat_cdf_persample['n_reads'], dat_cdf_persample[sn])

    dat_cdf_persample = pd.DataFrame(dat_cdf_persample)
    ymin = dat_cdf_persample.values[:,1:].min()
    ymax = dat_cdf_persample.values[:,1:].max()
    plt.plot([500, 500], [ymin * 0.8, ymax * 1.05], color='black')
    plt.ylabel('# Cells with >X reads')
    plt.xlabel('# Reads assigned to cell')

    plt.tight_layout()
    pdf.savefig()
    plt.figure()

    i_nr_closest = np.argmin(np.abs(dat_cdf_persample.n_reads-500))
    nr_closest = dat_cdf_persample.n_reads.values[i_nr_closest]

    plt.bar(dat_cdf_persample.columns[1:], dat_cdf_persample.values[i_nr_closest,1:])
    plt.xticks(rotation=90);
    plt.xlabel('Sample'); plt.ylabel('# Of Cells with >= %d reads' % nr_closest)
    plt.rc('axes', titlesize=10, labelsize=8)
    plt.tick_params(labelsize=6)

    plt.tight_layout()
    pdf.savefig()
    plt.figure()
   
    expected_per_60mil = 2500
    adjust_factor = float(library_total_reads) / 60000000
    good_covg = int(0.5 + expected_per_60mil * adjust_factor)


    i_nr_closest = np.argmin(np.abs(dat_cdf_persample.n_reads-good_covg))
    nr_closest = dat_cdf_persample.n_reads.values[i_nr_closest]
    
    print('EP60: %d, adjust: %f, final: %d, nr_c: %d' % (expected_per_60mil, adjust_factor, good_covg, nr_closest))

    plt.bar(dat_cdf_persample.columns[1:], dat_cdf_persample.values[i_nr_closest,1:])
    plt.xticks(rotation=90);
    plt.xlabel('Sample'); plt.ylabel('# Of Cells with >= %d Reads\n(Well-covered, adjusted for total reads)' % nr_closest)
    plt.rc('axes', titlesize=10, labelsize=8)
    plt.tick_params(labelsize=6)

    plt.tight_layout()
    pdf.savefig()
    plt.figure()

    ncell_samples = dat[np.logical_and(
                          dat.reads >= 10, 
                          dat.reads < nr_closest)].groupby('sample')['reads'].count().reset_index()
    plt.bar(ncell_samples['sample'], ncell_samples.reads)
    plt.xticks(rotation=90);
    plt.xlabel('Sample'); plt.ylabel('# Of Cells with < %d reads' % nr_closest)
    plt.rc('axes', titlesize=10, labelsize=8)
    plt.tick_params(labelsize=6)

    plt.tight_layout()
    pdf.savefig()



    dat = records[records.cell != '*']
    dat = dat[dat.assignment_valid == 'valid']

    plate_list = [1, 2, 3, 4]
    print(dat.well1_plate.value_counts())
    print(dat.well2_plate.value_counts())
    # get the maximum reads in a r/c for scale
    reads_max_w1 = dat.groupby(['well1_plate', 'well1_row', 'well1_col'])['reads'].sum().reset_index().reads.max()
    reads_max_w2 = dat.groupby(['well1_plate', 'well1_row', 'well1_col'])['reads'].sum().reset_index().reads.max()

    for plate in plate_list:
        hm_dat = dat[(dat.well1_plate == plate)].groupby(['well1_row', 'well1_col'])['reads'].sum().reset_index()
        hm_dat = hm_dat.pivot(index='well1_row', columns='well1_col', values='reads').fillna(0)
        plt.figure(figsize=(12,12))
        prop = 100 * hm_dat/hm_dat.values.sum()
        sbn.heatmap(hm_dat, annot=prop, vmin=0, vmax=reads_max_w1)
        title('Split 1, Plate %d, Total Reads' % plate)
        plt.tight_layout()
        pdf.savefig()
       
        hm_dat = dat[(dat.well2_plate == plate)].groupby(['well2_row', 'well2_col'])['reads'].sum().reset_index()
        hm_dat = hm_dat.pivot(index='well2_row', columns='well2_col', values='reads').fillna(0)
        plt.figure(figsize=(12,12))
        prop = 100 * hm_dat/hm_dat.values.sum()
        sbn.heatmap(hm_dat, annot=prop, vmin=0, vmax=reads_max_w2)
        title('Split 2, Plate %d, Total Reads' % plate)
        plt.tight_layout()
        pdf.savefig()

    dat.loc[:, 'well1'] = 96 * (dat.well1_plate - 1) + 12 * (dat.well1_row - 1) + (dat.well1_col - 1) + 1
    dat.loc[:, 'well2'] = 96 * (dat.well2_plate - 1) + 12 * (dat.well2_row - 1) + (dat.well2_col - 1) + 1
    # now the big-arse heatmaps
    hm_dat = dat.groupby(['well1', 'well2'])['reads'].sum().reset_index()
    hm_dat = hm_dat.pivot(index='well1', columns='well2', values='reads').fillna(0)
    plt.figure(figsize=(16,16))
    sbn.heatmap(hm_dat, annot=False, cmap='nipy_spectral')
    title('Total Reads: Split Wells')
    plt.tight_layout()
    pdf.savefig()

    hm_max = dat.groupby(['well1', 'well2'])['reads'].max().reset_index()
    hm_max = hm_max.pivot(index='well1', columns='well2', values='reads').fillna(0)
    plt.figure(figsize=(16,16))
    sbn.heatmap(hm_max, annot=False, cmap='nipy_spectral')
    title('Maximum Reads / Cell: Split Wells')
    plt.tight_layout()
    pdf.savefig()

    plt.figure(figsize=(16,16))
    sbn.heatmap(100 * (1 + hm_max)/(1 + hm_dat), annot=False, cmap='coolwarm')
    title('% of reads belonging to best cell: Split wells')
    plt.tight_layout()
    pdf.savefig()

    plt.figure(figsize=(16,16))
    hm_dat = dat[dat.reads >= 2500].groupby(['well1', 'well2'])['reads'].count().reset_index()
    if hm_dat.shape[0] >= 5:
        hm_dat = hm_dat.pivot(index='well1', columns='well2', values='reads').fillna(0)
        sbn.heatmap(hm_dat, annot=False, cmap='gist_heat')
        title('# Of samples achieiving >= 2,500 reads: Split wells')
        plt.tight_layout()
        pdf.savefig()

    pdf.close()


if __name__ == '__main__':
    main(get_args())
