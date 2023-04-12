# Make saturation curves given the .bam files in a pipeline output directory
import pysam
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from multiprocessing import Pool
from collections import Counter
from csv import DictReader
import os

READ_PROFILES = np.hstack([np.arange(0.0025, 1, 1/75), [1]])  # 75 points; plus [1]

from argparse import ArgumentParser

def get_args():
    parser = ArgumentParser()
    parser.add_argument('pipeline_dir', help='The pipeline output directory holding the aligned merged lysis bam files')
    parser.add_argument('library_digest', help='The library digest CSV file')
    parser.add_argument('sample_digest', help='The sample digest CSV file')
    parser.add_argument('--prefix', help='The prefix for the output file names', default='MolecularSaturation')
    parser.add_argument('--threads', help='The number of threads to use', default=12, type=int)
    parser.add_argument('--efficiency_factor', help='The efficiency factor to use', default='R2A')

    return parser.parse_args()


def lib_saturation(bampath, probs):
    hdl = pysam.AlignmentFile(bampath)
    tags = [set() for _ in probs]
    umi_counts = np.zeros(len(tags), dtype=int)
    decoded_umi_counts = np.zeros(len(tags), dtype=int)
    read_counts = np.zeros(len(tags), dtype=int)
    decoded_read_counts = np.zeros(len(tags), dtype=int)
    aligned_read_counts = np.zeros(len(tags), dtype=int)
    umi_cell_counts = [Counter() for _ in probs]
    pstart = 0
    k_clear = 0
    for read in hdl:
        if read.is_supplementary or read.is_secondary:
            continue
        # check if we can clear
        if read.reference_start != pstart:
            k_clear += 1
            if k_clear > 1000:  # referesh every kb
                for j, tagset in enumerate(tags):
                    umi_counts[j] += len(tagset)
                    decoded_umi_counts[j] += len({u for u in tagset if '*' not in u})
                    umi_cell_counts[j].update((k.split('@')[-1] for k in tagset))
                tags = [set() for _ in probs]
                k_clear = 0
            pstart = read.reference_start
        p = np.random.random()
        idx = np.where(probs > p)[0]
        mi = '@' + (str(read.get_tag('MI')) if read.has_tag('MI') else 'x') + '@'
        cb = read.get_tag('CB') if read.has_tag('CB') else '*'
        cb = '*' if cb == '' else cb
        loc = str(read.reference_start) if read.is_mapped else '*'
        tag = loc + mi + cb
        for j in idx:
            read_counts[j] += 1
            tags[j].add(tag)
        if cb != '*':
            for j in idx:
                decoded_read_counts[j] += 1
        if read.is_mapped and cb != '*':
            for j in idx:
                aligned_read_counts[j] += 1
            
    for j, tagset in enumerate(tags):
        umi_counts[j] += len(tagset)
        decoded_umi_counts[j] += len({u for u in tagset if '*' not in u})
        umi_cell_counts[j].update((k.split('@')[-1] for k in tagset))
    hdl.close()
    # want to count the UMI efficiency (foreground vs background)
    # typically we have ~5000 cells per sublibrary, and therefore
    # would expect 1/5000 UMI per cell. However the reads typically
    # follow some kind of Pareto distribution where even log counts
    # are right tailed. A 1/10000 rate of UMI (10e-4) could be a
    # plausible threshold; and to be conservative we use 1/20000
    # as there may be as few as 3000 cells in the sublibrary.
    background_umi = 0 * umi_counts
    for i, (tot_umi, cell_ctr) in enumerate(zip(umi_counts, umi_cell_counts)):
        thr = 1 + int(5e-5 * tot_umi)  # round to at least 1 for low-cvg libraries
        background_umi[i] = sum((v for v in cell_ctr.values() if v <= thr))
    return (read_counts, decoded_read_counts, 
            aligned_read_counts, umi_counts, 
            decoded_umi_counts, background_umi)


def lib_saturation_(args):
    return lib_saturation(*args)

def main(args):
    # want to map sublibraries to lysis IDs
    seq2lys = dict()
    seq2typ = dict()
    for record in DictReader(open(args.library_digest, 'rt')):
        seq2lys[record['sequence_id']] = record['lysis_id']
        seq2typ[record['sequence_id']] = record['library_type']
    assay2sample = dict()
    for record in DictReader(open(args.sample_digest, 'rt')):
        assay2sample[record['assay_info']] = record['sample_id']
    bam_files = [args.pipeline_dir + '/' + x for x in os.listdir(args.pipeline_dir) if x.endswith('.bam')]
    print(bam_files)
    sat_args = [(bam, READ_PROFILES) for bam in bam_files]
    pool = Pool(args.threads)
    coverages = pool.map(lib_saturation_, sat_args)
    pool.close()
    records = list()
    for bam, (nread, dread, aread, numi, dumi, bumi) in zip(bam_files, coverages):
        bam_name = bam.split('/')[-1]
        sublib = bam_name.split('__')[0]
        if '__UNK__' in bam_name:
            lysis = bam_name.split('__')[0].split('_')[0]
            type_ = 'rna' if 'Aligned' in bam_name else 'dna'
        else:
            lysis = seq2lys[sublib]
            type_ = seq2typ[sublib]
        sample = bam_name.split('__')[1]
        sample_name = assay2sample.get(sample, 'UNK')
        antibody = bam_name.split('__')[2]
        for frac, nr, dr, ar, nu, du, bu in zip(READ_PROFILES, nread, dread,
                                                aread, numi, dumi, bumi):
            records.append({
                'bam': bam_name,
                'sublibrary': sublib,
                'lysis': lysis,
                'library_type': type_,
                'sample': sample,
                'sample_name': sample_name,
                'antibody': antibody,
                'fraction': frac,
                'reads': nr,
                'decoded_reads': dr,
                'aligned_reads': ar,
                'umi': nu,
                'decoded_umi': du,
                'background_umi': bu,
                'molecular_efficiency': nu/nr,
                'efficiency': du/dr,
                'background_rate': bu/du
            })
    records = pd.DataFrame.from_records(records)
    print(records.head().to_string())
    records.to_csv(args.pipeline_dir + '/' + args.prefix + '.csv', index=False)

    postprocess_efficiency(records, args.pipeline_dir + '/' + args.prefix, USE_=args.efficiency_factor)


def postprocess_efficiency(erecords, outbase, USE_='R2A'):
    DCOLS = ['#2f4f4f',
             '#8b4513',
             '#228b22',
             '#00008b',
             '#ff0000',
             '#ffff00',
             '#00ff00',
             '#00ffff',
             '#ff00ff',
             '#6495ed',
             '#ff69b4',
             '#ffe4c4']
    def neglogitcdf(x, m, s):
        return 1/(1+np.exp(-(-x-m)/s))
    
    min_ = np.log10(100000)
    max_ = np.log10(3.5*10**9)
    reads_extrapolated = 10**np.arange(min_, max_, (max_-min_)/125)
    combined_info = None
    adjustments = dict()
    edata = erecords[erecords.sublibrary.map(lambda s: 'UNK' not in s)].copy()
    udata = erecords[erecords.sublibrary.map(lambda s: 'UNK' in s)].copy()
    
    for lys in edata.sublibrary.unique():
        ldat = edata[edata.sublibrary == lys]
        if USE_ == 'RAW':
            libdata = ldat.groupby('fraction')[['reads', 'umi']].apply(sum)
        elif USE_ == 'R2A':
            libdata = ldat.groupby('fraction')[['reads', 'decoded_umi']].apply(sum)
            libdata.columns = ['reads', 'umi']
        else:
            libdata = ldat.groupby('fraction')[['aligned_reads', 'decoded_umi']].apply(sum)
            libdata.columns = ['reads', 'umi']
        libdata['efficiency'] = libdata.umi/libdata.reads
        eadj = libdata.efficiency.max()
        adjustments[lys] = eadj
        # renormalized
        libdata['efficiency'] = libdata['efficiency']/eadj
        libdata['method'] = 'empirical'
        xfit = np.log(np.array([1] + libdata.reads.values.tolist()))
        yfit = np.array([1] + libdata.efficiency.values.tolist())
        epar, ecov = sp.optimize.curve_fit(neglogitcdf, xfit, yfit, method='lm', maxfev=10**5)
        epred = neglogitcdf(np.log(reads_extrapolated), *epar)
        extr_data = pd.DataFrame.from_dict({'reads': reads_extrapolated, 'umi': epred * reads_extrapolated * eadj, 'efficiency': epred})
        extr_data['method'] = 'extrapolated'
        mg_data = pd.concat([libdata.reset_index()[['reads', 'umi', 'efficiency', 'method']], extr_data[['reads', 'umi', 'efficiency', 'method']]], axis=0)
        mg_data['sublibrary'] = lys
        if combined_info is None:
            combined_info = mg_data
        else:
            combined_info = pd.concat([combined_info, mg_data])
    
    lysis_info = dict()
    for i, sublib in enumerate(combined_info.sublibrary.unique()):
        msub = combined_info[(combined_info.sublibrary == sublib) & (combined_info.method == 'extrapolated')].copy()
        msub.loc[:, 'umi_delta'] = np.hstack([[msub.umi.values[0]], msub.umi.values[1:] - msub.umi.values[:-1]])
        msub.loc[:, 'read_delta'] = np.hstack([[msub.reads.values[0]], msub.reads.values[1:] - msub.reads.values[:-1]])
        msub.loc[:, 'delta_ratio'] = msub.umi_delta / msub.read_delta
        target_idx = np.where(msub.delta_ratio.values < 0.05)[0]
        if len(target_idx) > 0:
            target_idx = target_idx[0]
        else:
            target_idx = msub.delta_ratio.values.shape[0] - 1
        lysis_id = edata[edata.sublibrary == sublib].lysis.values[0]
        library_type = edata[edata.sublibrary == sublib].library_type.values[0]
        if lysis_id not in lysis_info:
            lysis_info[lysis_id] = {'lysis_id': lysis_id}
        
        # pull in the rawdata
        esub = edata[(edata.sublibrary == sublib) & (edata.fraction == 1.)]
        usub = udata[(udata.lysis == lysis_id) & (udata.library_type == library_type) & (udata.fraction == 1.)]
        lysis_info[lysis_id][f'{library_type}_sequenced_reads'] = usub.reads.sum() + esub.reads.sum()
        lysis_info[lysis_id][f'{library_type}_reads'] = esub.reads.sum()
        lysis_info[lysis_id][f'{library_type}_barcode_efficiency'] = '%.1f%%' % (100 * esub.reads.sum()/(usub.reads.sum() + esub.reads.sum()))
        lysis_info[lysis_id][f'{library_type}_aligned_reads'] = esub.aligned_reads.sum()
        lysis_info[lysis_id][f'{library_type}_umi'] = esub.umi.sum()
        lysis_info[lysis_id][f'{library_type}_aligned_umi'] = esub.decoded_umi.sum()
        lysis_info[lysis_id][f'{library_type}_molecular_efficiency'] = '%.1f%%' % (100*lysis_info[lysis_id][f'{library_type}_umi']/lysis_info[lysis_id][f'{library_type}_reads'])
        lysis_info[lysis_id][f'{library_type}_efficiency'] = '%.1f%%' % (100 * lysis_info[lysis_id][f'{library_type}_aligned_umi']/lysis_info[lysis_id][f'{library_type}_aligned_reads'])
        lysis_info[lysis_id][f'{library_type}_effective_library_size'] = '%dM' % (int(msub.umi.values[target_idx]/10**6))
        lysis_info[lysis_id][f'{library_type}_reads_to_target'] = '%dM' % (int(msub.reads.values[target_idx]/10**6))
        lysis_info[lysis_id][f'{library_type}_efficiency_factor'] = '%.1f%%' % (100 * adjustments[sublib])
        if lysis_info[lysis_id][f'{library_type}_effective_library_size'] == '0M':
            lysis_info[lysis_id][f'{library_type}_effective_library_size'] = '<1M'
            

    lysis_records = pd.DataFrame(lysis_info).T.reset_index().iloc[:,1:]
    orig_ix = [x[4:] if x != 'lysis_id' else x for x in lysis_records.columns]
    lysis_records = lysis_records[sorted(lysis_records.columns, key = lambda x: orig_ix.index(x[4:]) if x != 'lysis_id' else 0)]
    lysis_records.to_csv(outbase + '_library_efficiency.csv', index=False)
    
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(outbase + '_efficiency_plots.pdf') as pdf:
        xlim = (combined_info.reads.min(), combined_info.reads.max())
        ylim = (combined_info.umi.min(), combined_info.umi.max())
        fig, ax = plt.subplots()
        for i, sublib in enumerate(combined_info.sublibrary.unique()):
            esub = combined_info[(combined_info.sublibrary == sublib) & (combined_info.method == 'extrapolated')]
            plt.plot(esub.reads/10**6, esub.umi/10**6, c=DCOLS[i], linestyle='dashed', label=f'{sublib}, proj')
            plt.legend()
            plt.xlabel('Decoded, Aligned Reads (M)')
            plt.ylabel('UMI (M)')
        pdf.savefig(fig)
        plt.xlim([0.1, 10**3])
        newy=combined_info[combined_info.reads/10**6 < 10**3].umi.max()*1.05/10**6
        plt.ylim([ylim[0]/10**6, newy])
        pdf.savefig(fig)
        plt.xlim([0.1, 2.5*10**2])
        newy=combined_info[combined_info.reads/10**6 < 2.5*10**2].umi.max()*1.05/10**6
        plt.ylim([ylim[0]/10**6, newy])
        pdf.savefig(fig)
        plt.close()
        
        fig, ax = plt.subplots()
        for i, sublib in enumerate(combined_info.sublibrary.unique()):
            msub = combined_info[(combined_info.sublibrary == sublib) & (combined_info.method == 'empirical')]
            plt.plot(msub.reads/10**6, msub.efficiency * adjustments[sublib], c=DCOLS[i], linestyle='solid', label=f'{sublib}, raw')
            esub = combined_info[(combined_info.sublibrary == sublib) & (combined_info.method == 'extrapolated')]
            plt.plot(esub.reads/10**6, esub.efficiency * adjustments[sublib], c=DCOLS[i], linestyle='dashed', label=f'{sublib}, proj')
        plt.legend()
        plt.xlabel('Decoded, Aligned Reads (M)')
        plt.ylabel('Efficiency (Analyzable UMI/Rd)')
        xlim=[combined_info[combined_info.method == 'empirical'].reads.min()/10**6, combined_info[combined_info.method == 'empirical'].reads.max()/10**5]
        plt.xlim(xlim)
        ax.set_xscale('log')
        pdf.savefig(fig)
        plt.close()
    


if __name__ == '__main__':
    main(get_args())
