import pandas as pd
import seaborn as sbn
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from argparse import ArgumentParser
from csv import DictReader
import warnings
warnings.filterwarnings('ignore')

COLOR_MAP = ['orange', 'purple', 'blue', 'red', 'black']
MARKER_MAP = ['.', '+', '1', 's']

ER_QUANTILES = (5, 95)
MIN_READS_QUANTILE = 0.01

parser = ArgumentParser()
parser.add_argument('cell_qc', help='The cell qc file')
parser.add_argument('sample_qc', help='The sample qc file')
parser.add_argument('out_pdf', help='The output pdf')

args = parser.parse_args()

print('reading...')
chip_cell_qc = pd.read_csv(args.cell_qc, sep='\t')
chip_sample_qc = pd.read_csv(args.sample_qc, sep='\t')


print('cell size: %s, sample_size: %s' % (repr(chip_cell_qc.shape), repr(chip_sample_qc.shape)))

def unmelt(df, cond_list, value_column, ix_columns):
    um_df = None
    for cols, vals in cond_list:
        print(' .. .. unmelting %s : %s key by %s' % (repr(cols), repr(vals), repr(ix_columns)))
        df_sub = df.copy()
        keyname = '_'.join(vals)
        for col, val in zip(cols, vals):
            df_sub = df_sub[df_sub.loc[:, col] == val]
        df_sub.index = df_sub.loc[:, ix_columns].agg(':'.join, axis=1)
        srs = df_sub.loc[:, value_column]
        if um_df is None:
            um_df = df_sub.copy()
        srsd = {k: v for k, v in zip(srs.index, srs.tolist())}
        um_df.loc[:, keyname] = np.array([srsd[x] for x in um_df.index])
    return um_df


with PdfPages(args.out_pdf) as pdf:
    print('.. computing read statistics ..')
    read_counts = unmelt(chip_cell_qc, [(('filter', 'target', 'type'), ('retained', 'off_target', 'reads')),
                                        (('filter', 'target', 'type'), ('retained', 'in_peak', 'reads'))],
                       'count', ('library', 'assay', 'cell'))
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
    
    cell_stats = unmelt(chip_cell_qc, [
        (('filter', 'target', 'type'), ('retained', 'in_peak', 'umi')),
        (('filter', 'target', 'type'), ('retained', 'off_target', 'umi')),
        (('filter', 'target', 'type'), ('retained', 'in_enhancer', 'umi')),
        (('filter', 'target', 'type'), ('retained', 'in_promoter', 'umi')),
        (('filter', 'target', 'type'), ('retained', 'in_genebody', 'umi')),
        (('filter', 'target', 'type'), ('retained', 'in_peak', 'reads')),
        (('filter', 'target', 'type'), ('retained', 'off_target', 'reads'))
        ], 'count', ('library', 'assay', 'cell'))
    cell_stats = cell_stats[cell_stats.index.isin(good_cells)]
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
    plt.gca().set_rasterized(True)
    #plt.savefig(pdf, format='pdf', dpi=250)
    plt.close()
    plt.figure()
    print(' .. .. umi vs reads/umi')
    plt.scatter(cell_stats.total_umi, cell_stats.reads_per_umi, marker='.')
    plt.xlabel('Total UMI'); plt.ylabel('Reads per UMI'); plt.title('Reads per UMI by Total UMI (per cell)\nAll libraries')
    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    plt.gca().set_rasterized(True)
    #plt.savefig(pdf, format='pdf', dpi=250)
    plt.close()
    plt.figure()
    print(' .. .. umi vs reads/umi (colored)')
    # pass 1: compute the mean reads / umi
    rpu = cell_stats[cell_stats.total_umi >= 100].groupby('library')['reads_per_umi'].mean().reset_index()
    new_names = {lib: '%s [RPU=%.2f]' % (lib, rpu) for lib, rpu in zip(rpu.library, rpu.reads_per_umi)}
    libs_sorted = rpu.library[np.argsort(rpu.reads_per_umi)]
    umi_min=cell_stats.total_umi.min()
    umi_max=cell_stats.total_umi.max()
    for lib in libs_sorted:
        print(lib)
        csub=cell_stats[(cell_stats.total_umi > 100) & (cell_stats.library == lib)].copy()
        print(csub)
        plt.scatter(csub.total_umi, 1./csub.reads_per_umi, marker='.', label=new_names[lib])
        if csub.shape[0] > 0:
            mean_rpu = np.mean(1./csub.reads_per_umi)
        plt.plot((umi_min, umi_max), (mean_rpu, mean_rpu))
    plt.xlabel('Total UMI in cell'); plt.ylabel('UMI per read per cell'); plt.title('Reads per UMI by Total UMI (per cell)\nAll libraries')
    #plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    plt.legend(bbox_to_anchor=(1,1.00))
    plt.tight_layout()
    plt.gca().set_rasterized(True)
    #plt.savefig(pdf, format='pdf', dpi=250)
    plt.close()
    plt.figure(figsize=(12,8))
    print(' .. .. reads per umi : assay boxplot')
    print(cell_stats.head())
    try:
      sbn.boxplot(x='assay', y='reads_per_umi', hue='library', data=cell_stats[cell_stats.total_umi >= 100])
      plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True)
      #plt.savefig(pdf,format='pdf',dpi=250)
      plt.close()
    except:
      pass
    
    plt.figure(figsize=(8,6))
    cell_stats.loc[:,'umi_per_read'] = 1/cell_stats.reads_per_umi
    cell_stats.loc[:,'total_reads_log'] = np.log10(cell_stats.total_reads)
    cell_stats.loc[:, 'total_umi_log'] = np.log10(cell_stats.total_umi)
    sbn.lmplot(x='total_reads_log', y='umi_per_read', hue='library',
               data=cell_stats)
    plt.xlabel('Total Reads (log10)'); plt.ylabel('UMI per Read')
    plt.gca().set_rasterized(True)
    #plt.savefig(pdf, format='pdf', dpi=250)
    plt.close()
    
    plt.figure(figsize=(12,8))
    sbn.lmplot(x='total_umi_log', y='umi_per_read', hue='library',
               data=cell_stats)
    plt.xlabel('Total UMI (log10)'); plt.ylabel('UMI per Read')
    plt.gca().set_rasterized(True)
    #plt.savefig(pdf, format='pdf', dpi=250)
    plt.close()
    
    plt.figure(figsize=(8,6))
    print('.. .. total umi per cell, violin')
    sbn.violinplot(x='assay', y='total_umi', hue='library', data=cell_stats)
    plt.legend(bbox_to_anchor=(1.04,1))
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True)
    #plt.savefig(pdf,format='pdf',dpi=250)
    plt.close()
    
    plt.figure(figsize=(12,8))
    print('.. .. total umi per cell, boxplot')
    cell_stats.loc[:, 'total_umi_log10'] = np.log10(cell_stats.total_umi)
    sbn.boxplot(x='assay', y='total_umi_log10', hue='library', data=cell_stats)
    plt.legend(bbox_to_anchor=(1, 1.04))
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True)
    #plt.savefig(pdf,format='pdf',dpi=250)
    plt.close()
    
    print('.. .. swarmplot by library')
    for lib in cell_stats.library.unique():
        dat = cell_stats[cell_stats.library == lib]
        print('.. .. .. %s : %d' % (lib, dat.shape[0]))
        if dat.shape[0] > 15000:
            print('.. .. .. More than 15k cells for library %s -- subsetting' % lib)
            dat.loc[:, 'srt'] = dat.total_umi_log10.values + 0.01 * np.random.uniform(size=(dat.shape[0],))
            dat = dat[np.argsort(np.argsort(-dat.srt)) < 15000]
            print('.. .. .. done: %d' % dat.shape[0])
        plt.figure(figsize=(12,8))
        sbn.swarmplot(x='assay', y='total_umi_log10', data=dat)
        plt.title('Cell UMI for library %s' % lib)
        plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True)
        #plt.savefig(pdf,format='pdf',dpi=250)
        plt.close()
        plt.figure()
    
    cell_stats.loc[:, 'enhancer_rate'] = cell_stats.retained_in_enhancer_umi / cell_stats.total_umi
    cell_stats.loc[:, 'promoter_rate'] = cell_stats.retained_in_promoter_umi / cell_stats.total_umi
    cell_stats.loc[:, 'retained_in_CRE_umi'] = cell_stats.retained_in_enhancer_umi + cell_stats.retained_in_promoter_umi
    cell_stats.loc[:, 'CRE_rate'] = cell_stats.retained_in_CRE_umi / cell_stats.total_umi
    cell_stats.loc[:, 'genic_rate'] = cell_stats.retained_in_genebody_umi / cell_stats.total_umi
 
    cell_stats.to_csv(args.out_pdf[:-4] + '_cell_stats.csv')

    plt.figure(figsize=(8,8))
    for ab in cell_stats.antibody.unique():
        print('.. ab enhancer/promoter rates ..' + str(ab))
        dsub = cell_stats[cell_stats.antibody == ab]
        plt.scatter(dsub.enhancer_rate, dsub.promoter_rate, label=ab, marker='.')
    plt.xlabel('Enhancer Rate (ENCODE pELS/dELS UMI / total)');
    plt.ylabel('Promoter Rate (ENCODE PLS UMI / total)');
    plt.title('Enhancer / Promoter rates (all libraries)');
    plt.legend();
    plt.gca().set_rasterized(True)
    #plt.savefig(pdf, format='pdf', dpi=250)
    plt.close()
   
    print('.. library swarms ..') 
    for lib in cell_stats.library.unique():
        dat = cell_stats[cell_stats.library == lib]
        if dat.shape[0] > 15000:
            print('.. .. .. More than 15k cells for library %s -- subsetting' % lib)
            dat.loc[:, 'srt'] = dat.total_umi_log10.values + 0.01 * np.random.uniform(size=(dat.shape[0],))
            dat = dat[np.argsort(np.argsort(-dat.srt)) < 15000]
        plt.figure(figsize=(12,8))
        sbn.swarmplot(x='assay', y='CRE_rate', hue='antibody', data=dat)
        plt.title('Rate of UMI in ENCODE-defined CRE\n(Library %s)' % lib)
        plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True)
        #plt.savefig(pdf,format='pdf',dpi=250)
        plt.close()
        plt.figure(figsize=(12,8))
        sbn.swarmplot(x='assay', y='genic_rate', hue='antibody', data=dat)
        plt.title('Rate of UMI in genic regions\n(Library %s)' % lib)
        plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True)
        #plt.savefig(pdf,format='pdf',dpi=250)
        plt.close()
    
   
    print(' .. antibody scatters ..') 
    # make colors for the antibodies
    plotted=False
    plt.figure()
    for i, library in enumerate(cell_stats.library.unique()):
        if (i % len(MARKER_MAP)) == 0:
            if plotted:
                plt.legend(bbox_to_anchor=(1.04,1.04))
                #pdf.savefig()
                plotted = False
        dat = cell_stats[cell_stats.library == library]
        for j, ab in enumerate(dat.antibody.unique()):
            dat2 = dat[dat.antibody == ab]
            if not plotted:
                plt.figure(figsize=(12,8))
            plt.scatter(dat2.total_umi_log10, dat2.CRE_rate, label=str(ab) + ':' + str(library),
                        color=COLOR_MAP[j % len(COLOR_MAP)],
                        marker=MARKER_MAP[i % len(MARKER_MAP)])
            plt.xlabel('Total UMI (log10)')
            plt.ylabel('CRE rate')
            plotted = True
    
    if plotted:
        plt.legend(bbox_to_anchor=(1.04,1.04))
        #pdf.savefig()
        plt.close()
        
    print('.. computing sample stats ..') 
    sample_dat = unmelt(chip_sample_qc,
                       [(('feature', 'metric'), ('HQ_umi', 'total')),
                        (('feature', 'metric'), ('HQ_umi', 'median')),
                        (('feature', 'metric'), ('umi_CRE_rate', 'median')),
                        (('feature', 'metric'), ('umi_genic_rate', 'median')),
                        (('feature', 'metric'), ('enhancer_promoter_ratio', 'total'))],
                       'value', ('sample_id',))
   
    print('..plotting sample stats..') 
    plt.figure()
    sbn.barplot(x='sample', y='HQ_umi_total', hue='antibody', data=sample_dat)
    plt.legend(bbox_to_anchor=(1.04,1.04))
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
    
    plt.figure(figsize=(12,8))
    sbn.barplot(x='sample', y='HQ_umi_median', hue='antibody', data=sample_dat)
    plt.legend(bbox_to_anchor=(1.04,1.04))
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
    
    plt.figure(figsize=(12,8))
    sbn.barplot(x='sample', y='umi_CRE_rate_median', hue='antibody', data=sample_dat)
    plt.legend(bbox_to_anchor=(1.04,1.04))
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
    
    plt.figure(figsize=(12,8))
    sbn.barplot(x='sample', y='umi_genic_rate_median', hue='antibody', data=sample_dat)
    plt.legend(bbox_to_anchor=(1.04,1.04))
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)
    
    plt.figure(figsize=(12,8))
    sbn.barplot(x='sample', y='enhancer_promoter_ratio_total', hue='antibody', data=sample_dat)
    plt.legend(bbox_to_anchor=(1.04,1.04));
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True);plt.savefig(pdf,format='pdf',dpi=250)

    sample_dat.to_csv(args.out_pdf[:-4] + '_sample_stats.csv')

    # lastly: extract the unmapped information per library and sample
    reader = DictReader(open(args.cell_qc, 'rt'), delimiter='\t')
    unm_counts = dict()
    for entry in reader:
        if entry['target'] in {'in_peak', 'off_target'} and entry['type'] == 'reads' and entry['filter'] == 'retained':
            if entry['sample'] == '*' or entry['cell'] == '*':
                lib = entry['library']
                if lib not in unm_counts:
                    unm_counts[lib] = dict()
                if entry['assay'] not in unm_counts[lib]:
                    unm_counts[lib][entry['assay']] = 0
                unm_counts[lib][entry['assay']] += int(entry['count'])
    unm_dat = list()
    for lib, cdc in unm_counts.items():
        for assay, cnt in cdc.items():
            unm_dat.append({'library': lib, 'assay': assay, 'reads': cnt})
    unm_dat = pd.DataFrame(unm_dat)

    plt.figure(figsize=(12,8))
    gp = sbn.barplot(x='library', y='reads', hue='assay', data=unm_dat)
    gp.set_yscale('log')
    plt.title('Unmatched reads')
    plt.xticks(rotation=90);plt.tight_layout();plt.gca().set_rasterized(True)
    #plt.savefig(pdf,format='pdf',dpi=250)
    #pdf.savefig()
             
        
