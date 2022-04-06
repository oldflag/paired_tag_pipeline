from argparse import ArgumentParser
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import seaborn as sbn
import numpy as np
import gzip

def get_args():
    parser = ArgumentParser()
    parser.add_argument('barcode_csv', help='The barcode csv.gz file')
    parser.add_argument('output_pdf', help='The output qc pdf')

    return parser.parse_args()


def main(args):
    hdl = gzip.open(args.barcode_csv, 'rt')
    counts = dict()
    for line in hdl:
        vals = line.strip().split(',')[1].split(':')
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

    r_unassign = records[records.sample_cell.map(lambda x: '*' in x)].reads.sum() / records.reads.sum()
    str_A = 'assigned (%.1f%%)' % (100 * (1-r_unassign)) 
    str_U = 'unassigned (%.1f%%)' % (100 * r_unassign)
    records.loc[:, 'assignment'] = records.sample_cell.map(lambda x: str_U if '*' in x else str_A)
    acounts = records.groupby('assignment')['reads'].sum().reset_index()
    pdf = PdfPages(args.output_pdf)
    plt.bar(acounts.assignment, acounts.reads)
    plt.ticklabel_format(axis='y', style='plain')

    plt.tight_layout()
    pdf.savefig()
    plt.figure()

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

    pdf.close()


if __name__ == '__main__':
    main(get_args())
