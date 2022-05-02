import pandas as pd
import seaborn as sbn
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from argparse import ArgumentParser
import warnings
warnings.filterwarnings('ignore')

COLOR_MAP = ['orange', 'purple', 'blue', 'red', 'black']
MARKER_MAP = ['.', '+', '1', 's']

ER_QUANTILES = (10, 99)
MIN_READS_QUANTILE = 0.05

parser = ArgumentParser()
parser.add_argument('cell_qc', help='The cell qc file')
parser.add_argument('sample_qc', help='The sample qc file')
parser.add_argument('out_pdf', help='The output pdf')

args = parser.parse_args()

print('reading...')
chip_cell_qc = pd.read_csv(args.cell_qc, sep='\t')
chip_sample_qc = pd.read_csv(args.sample_qc, sep='\t')


print('cell size: %s, sample_size: %s' % (repr(chip_cell_qc.shape), repr(chip_sample_qc.shape)))

def unmelt(df, cond_list, value_column, ix_column):
    um_df = None
    for cols, vals in cond_list:
        df_sub = df
        keyname = '_'.join(vals)
        for col, val in zip(cols, vals):
            df_sub = df_sub[df_sub.loc[:, col] == val]
        df_sub.index = df_sub[ix_column]
        srs = df_sub.loc[:, value_column]
        if um_df is None:
            um_df = df_sub.copy()
        um_df.loc[:, keyname] = srs[df_sub.index]
    return um_df


with PdfPages(args.out_pdf) as pdf:
    print('.. computing read statistics ..')
    read_counts = unmelt(chip_cell_qc, [(('filter', 'target', 'type'), ('retained', 'off_target', 'reads')),
                                        (('filter', 'target', 'type'), ('retained', 'in_peak', 'reads'))],
                       'count', 'cell')
    read_counts['total_reads'] = read_counts.retained_off_target_reads + read_counts.retained_in_peak_reads
    
    p_il, x0 = np.percentile, read_counts.total_reads[read_counts.total_reads > 0]
    l0 = p_il(x0, ER_QUANTILES[0])
    x0 = x0[x0 > l0]
    l1 = p_il(x0, ER_QUANTILES[1])
    EFFECTIVE_RANGE = (l0, l1)
    MIN_READS_CELL = EFFECTIVE_RANGE[0] + MIN_READS_QUANTILE * (EFFECTIVE_RANGE[1] - EFFECTIVE_RANGE[0])
    MIN_READS_CELL
    
    print('Effective range: %d - %d; Minimum reads: %.1f' % (EFFECTIVE_RANGE[0],
                                                             EFFECTIVE_RANGE[1],
                                                             MIN_READS_CELL))
    
    good_cells = read_counts[read_counts.total_reads > MIN_READS_CELL].index
    if good_cells.shape[0] > 500000:
        print('.. more than 500K cells; downsampling by reads to 500K for plotting ..')
        read_counts.loc[:,'rnd'] = np.random.uniform(size=(read_counts.shape[0],))
        read_counts.loc[:,'srt'] = read_counts.total_reads + read_counts.rnd
        good_cells = read_counts[np.argsort(np.argsort(-read_counts.srt)) <= 500000].index
    print(' .. computing cell stats from size of %s ..' % repr(good_cells.shape))
    
    cell_qc_good = chip_cell_qc[chip_cell_qc.cell.isin(good_cells)].copy()
    cell_qc_good
    
    cell_stats = unmelt(cell_qc_good, [
        (('filter', 'target', 'type'), ('retained', 'in_peak', 'umi')),
        (('filter', 'target', 'type'), ('retained', 'off_target', 'umi')),
        (('filter', 'target', 'type'), ('retained', 'in_enhancer', 'umi')),
        (('filter', 'target', 'type'), ('retained', 'in_promoter', 'umi')),
        (('filter', 'target', 'type'), ('retained', 'in_genebody', 'umi')),
        (('filter', 'target', 'type'), ('retained', 'in_peak', 'reads')),
        (('filter', 'target', 'type'), ('retained', 'off_target', 'reads'))
        ], 'count', 'cell')
    cell_stats.loc[:, 'total_reads'] = cell_stats.retained_in_peak_reads + cell_stats.retained_off_target_reads
    cell_stats.loc[:,'total_umi'] = cell_stats.retained_in_peak_umi + cell_stats.retained_off_target_umi
    cell_stats.loc[:, 'reads_per_umi'] = cell_stats.total_reads/cell_stats.total_umi
    print(' .. cell statistic plots ..')
    plt.figure()
    print(' .. .. reads vs umi')
    plt.scatter(cell_stats.total_reads, cell_stats.total_umi, marker='.')
    plt.xlabel('Total Reads'); plt.ylabel('Total UMI'); plt.title('Reads and UMI (per cell)\nAll libraries')
    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    pdf.savefig()
    plt.figure(figsize=(12,8))
    print(' .. .. reads per umi : sample boxplot')
    sbn.boxplot(x='sample', y='reads_per_umi', hue='library', data=cell_stats)
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
    pdf.savefig()
    plt.figure(figsize=(12,8))
    print('.. .. total umi per cell, violin')
    sbn.violinplot(x='sample', y='total_umi', hue='library', data=cell_stats)
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
    plt.legend(bbox_to_anchor=(1.04,1))
    cell_stats.loc[:, 'total_umi_log10'] = np.log10(cell_stats.total_umi)
    pdf.savefig()
    plt.figure(figsize=(12,8))
    print('.. .. total umi per cell, boxplot')
    sbn.boxplot(x='sample', y='total_umi_log10', hue='library', data=cell_stats)
    plt.legend(bbox_to_anchor=(1, 1.04))
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
    pdf.savefig()
    print('.. .. swarmplot by library')
    for lib in cell_stats.library.unique():
        dat = cell_stats[cell_stats.library == lib]
        print('.. .. .. %s : %d' % (lib, dat.shape[0]))
        if dat.shape[0] > 50000:
            print('.. .. .. More than 15k cells for library %s -- subsetting' % lib)
            dat.loc[:, 'srt'] = dat.total_umi_log10.values + 0.01 * np.random.uniform(size=(dat.shape[0],))
            dat = dat[np.argsort(np.argsort(-dat.srt)) < 15000]
            print('.. .. .. done: %d' % dat.shape[0])
        plt.figure(figsize=(12,8))
        sbn.swarmplot(x='sample', y='total_umi_log10', data=dat)
        plt.title('Cell UMI for library %s' % lib)
        plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
        pdf.savefig()
    
    cell_stats.loc[:, 'enhancer_rate'] = cell_stats.retained_in_enhancer_umi / cell_stats.total_umi
    cell_stats.loc[:, 'promoter_rate'] = cell_stats.retained_in_promoter_umi / cell_stats.total_umi
    cell_stats.loc[:, 'retained_in_CRE_umi'] = cell_stats.retained_in_enhancer_umi + cell_stats.retained_in_promoter_umi
    cell_stats.loc[:, 'CRE_rate'] = cell_stats.retained_in_CRE_umi / cell_stats.total_umi
    cell_stats.loc[:, 'genic_rate'] = cell_stats.retained_in_genebody_umi / cell_stats.total_umi
    cell_stats.loc[:, 'antibody'] = cell_stats.sample_id.map(lambda x: x.split('__')[2])
   
    print('.. ab enhancer/promoter rates ..') 
    for ab in cell_stats.antibody.unique():
        dsub = cell_stats[cell_stats.antibody == ab]
        plt.figure(figsize=(8,8))
        plt.scatter(dsub.enhancer_rate, dsub.promoter_rate, label=ab, marker='.')
    plt.xlabel('Enhancer Rate (ENCODE pELS/dELS UMI / total)');
    plt.ylabel('Promoter Rate (ENCODE PLS UMI / total)');
    plt.title('Enhancer / Promoter rates (all libraries)');
    plt.legend();
    pdf.savefig()
   
    print('.. library swarms ..') 
    for lib in cell_stats.library.unique():
        dat = cell_stats[cell_stats.library == lib]
        if dat.shape[0] > 15000:
            print('.. .. .. More than 15k cells for library %s -- subsetting' % lib)
            dat.loc[:, 'srt'] = dat.total_umi_log10.values + 0.01 * np.random.uniform(size=(dat.shape[0],))
            dat = dat[np.argsort(np.argsort(-dat.srt)) < 15000]
        plt.figure(figsize=(12,8))
        sbn.swarmplot(x='sample', y='CRE_rate', hue='antibody', data=dat)
        plt.title('Rate of UMI in ENCODE-defined CRE\n(Library %s)' % lib)
        plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
        pdf.savefig()
        plt.figure(figsize=(12,8))
        sbn.swarmplot(x='sample', y='genic_rate', hue='antibody', data=dat)
        plt.title('Rate of UMI in genic regions\n(Library %s)' % lib)
        plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
        pdf.savefig()
    
   
    print(' .. antibody scatters ..') 
    # make colors for the antibodies
    plotted=False
    for i, library in enumerate(cell_stats.library.unique()):
        if (i % len(MARKER_MAP)) == 0:
            if plotted:
                plt.legend(bbox_to_anchor=(1.04,1.04))
                pdf.savefig()
                plotted = False
        dat = cell_stats[cell_stats.library == library]
        for j, ab in enumerate(dat.antibody.unique()):
            dat2 = dat[dat.antibody == ab]
            if not plotted:
                plt.figure(figsize=(12,8))
            plt.scatter(dat2.total_umi_log10, dat2.CRE_rate, label=ab + ':' + library,
                        color=COLOR_MAP[j % len(COLOR_MAP)],
                        marker=MARKER_MAP[i % len(MARKER_MAP)])
            plt.xlabel('Total UMI (log10)')
            plt.ylabel('CRE rate')
            plotted = True
    
    if plotted:
        plt.legend(bbox_to_anchor=(1.04,1.04))
        pdf.savefig()
        
    print('.. computing sample stats ..') 
    sample_dat = unmelt(chip_sample_qc,
                       [(('feature', 'metric'), ('HQ_umi', 'total')),
                        (('feature', 'metric'), ('HQ_umi', 'median')),
                        (('feature', 'metric'), ('umi_CRE_rate', 'median')),
                        (('feature', 'metric'), ('umi_genic_rate', 'median')),
                        (('feature', 'metric'), ('enhancer_promoter_ratio', 'total'))],
                       'value', 'sample_id')
   
    print('..plotting sample stats..') 
    plt.figure()
    sbn.barplot(x='sample', y='HQ_umi_total', hue='library', data=sample_dat)
    plt.legend(bbox_to_anchor=(1.04,1.04))
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
    pdf.savefig()
    plt.figure(figsize=(12,8))
    sbn.barplot(x='sample', y='HQ_umi_median', hue='library', data=sample_dat)
    plt.legend(bbox_to_anchor=(1.04,1.04))
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
    pdf.savefig()
    plt.figure(figsize=(12,8))
    sbn.barplot(x='sample', y='umi_CRE_rate_median', hue='library', data=sample_dat)
    plt.legend(bbox_to_anchor=(1.04,1.04))
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
    pdf.savefig()
    plt.figure(figsize=(12,8))
    sbn.barplot(x='sample', y='umi_genic_rate_median', hue='library', data=sample_dat)
    plt.legend(bbox_to_anchor=(1.04,1.04))
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
    pdf.savefig()
    plt.figure(figsize=(12,8))
    sbn.barplot(x='sample', y='enhancer_promoter_ratio_total', hue='library', data=sample_dat)
    plt.legend(bbox_to_anchor=(1.04,1.04));
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
    pdf.savefig()
