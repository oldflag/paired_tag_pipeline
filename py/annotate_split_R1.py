"""
This CLT annotates a Paired-Tag R1 fastq.gz with extracted barcodes from R2, and additionally
splits the file into mulitple output fastQs for use with parallel execution of the pipeline.

"""
from argparse import ArgumentParser
from utils import read_fasta, read_fastq, xopen, dxopen
from itertools import product, islice
import sys


def grouped_product(*iters, n):
    product_iter = product(*iters)
    entry = list(islice(product_iter, n))
    while entry:
           yield entry
           entry = list(islice(product_iter, n))


def setadd(st, x):
    st.add(x)
    return st


def get_args():
    parser = ArgumentParser('annotate_split_R1')
    parser.add_argument('R1_fastq', help='The Read 1 paired-tag fastQ, containing genomic sequence')
    parser.add_argument('bc_csv', help='The barcode CSV parsed from Read 2')
    parser.add_argument('well_bc', help='The well barcode file (fasta)')
    parser.add_argument('sample_bc', help='The sample barcode file (fasta)')
    parser.add_argument('--cells_per_fastq', help='The target # cells / fastQ for output', default=10000, type=int)
    parser.add_argument('--estimated_cells', help='An estimate of the total # of cells', type=int, default=None)
    parser.add_argument('--outdir', help='The output directory', default=None)

    return parser.parse_args()



def estimate_via_recapture(bcfile, n):
    """
    Use capture/recapture methodology to provide a coarse estimate
    of the number of cells per sample, number of samples, and total
    cells.

    Inputs
    --------
    :bcfile: The barcode file
    :n:      The number of (matching) reads to use in the estimate


    Outputs
    ---------
    number_of_cells, number_of_samples, cells_per_sample

    """
    nmatch = 0
    per_sample_umi = dict()
    # capture + mark
    with xopen(bcfile, 'rt') as handle:
        for entry in handle:
            bcinfo = entry.split(',')[1]
            if '*' not in bcinfo:
                codes = bcinfo.split(':')
                # umi:bc1:bc2:sm:re
                sample, cell, umi = codes[3], codes[1] + codes[2], codes[0]
                if sample not in per_sample_umi:
                    per_sample_umi[sample] = dict()
                per_sample_umi[sample][cell] = setadd(per_sample_umi[sample].get(cell, set()), codes[0])
                nmatch += 1
            if nmatch >= n/2:
                break
        # recapture
        per_sample_new_cells, per_sample_repeats = dict(), dict()
        nmatch = 0
        for entry in handle:
            bcinfo = entry.split(',')[1]
            if '*' not in bcinfo:
                codes = bcinfo.split(':')
                sample, cell, umi = codes[3], codes[1] + codes[2], codes[0]
                if sample not in per_sample_umi:
                    raise ValueError('Cell number estimation failed -- too few reads used')
                if cell not in per_sample_umi[sample]:
                    per_sample_new_cells[sample] = per_sample_new_cells.get(sample, 0) + 1
                    nmatch += 1
                else:
                    per_sample_repeats[sample] = per_sample_repeats.get(sample, 0) + \
                                                 (umi not in per_sample_umi[sample][cell])
                    nmatch += 1
            if nmatch >= n/2:
                break

    n_samples = len(per_sample_new_cells)
    cps = {sm: int(len(per_sample_umi[sm])*(per_sample_new_cells[sm] + per_sample_repeats[sm])/per_sample_repeats[sm])
            for sm in per_sample_umi}
    total_cells = sum(cps.values())
    return total_cells, n_samples, cps


def get_group_map(pool_fasta, sample_fasta, n_groups, out_base):
    """
    Produce a dictionary of output handles for groups of cell ids

    Input
    ---------
    :pool_fasta: Fasta file of pool barcodes
    :sample_fasta: Fasta file of sample barcodes
    :n_groups: Number of groups -- how many to divide the space

    Output
    ---------
    Map of {bc1:bc2:sm -> handle}, handle to (nomatch) reads

    """
    bc1_f, bc2_f = read_fasta(pool_fasta), read_fasta(pool_fasta)
    sm_f = read_fasta(sample_fasta)
    out_files = [f'{out_base}_{1+grp}.fq.gz' for grp in range(n_groups)]
    out_handles = [dxopen(u, 'wt') for u in out_files]
    handles = dict()
    for group_list in grouped_product(sm_f, bc1_f, bc2_f, n=n_groups):
        for i, (sm, bc1, bc2) in enumerate(group_list):
            handles[bc1.name + bc2.name + sm.name] = out_handles[i]

    sys.stderr.write(f'Writing to {len(out_handles) + 2} fastQ files\n') 

    return handles, dxopen('%s_unlinked.fq.gz' % out_base, 'wt'), dxopen('%s_sample_unlinked.fq.gz' % out_base, 'wt')
   

def get_base(fn, pfx=None):
    """
    Get the base of a fastq file by stripping:
      .gz .fq .fastq. .FQ .fastQ _R1 _1

    and optionally adding a prefix (dir)
    """
    pfx = '' if pfx is None else pfx.strip('/') + '/'  # ensure just one /
    if '_R1' in fn:
        bs = fn.split('/')[-1].split('_R1')[0]
    if '_1.f' in fn:
        bs = fn.split('/')[-1].split('_1.f')[0]
    else:
        bs = fn.strip('.gz').strip('.FQ').strip('.fastq').strip('.fq')  

    return pfx + bs
               

def main(args):
    if args.estimated_cells is None:
        ncell, _, _ = estimate_via_recapture(args.bc_csv, n=500000)  # use 500K reads to estimate # of cells
    else:
        ncell = args.estimated_cells


    n_groups = 1 + int(ncell / args.cells_per_fastq)
    sys.stderr.write(f'Using {ncell} estimated total cells and {n_groups} groups\n')

    base = get_base(args.R1_fastq, args.outdir)

    group_map, reject_handle, sample_only_handle = get_group_map(args.well_bc, args.sample_bc, n_groups, base)

    fq_in = read_fastq(args.R1_fastq)
    with xopen(args.bc_csv, 'rt') as bcinput:
        bc_in = (x.strip().split(',') for x in bcinput)

        for fq, bc in zip(fq_in, bc_in):
            # make sure that the names really do match
            assert fq.name.split()[0] == bc[0].split()[0], 'FastQ records and barcodes out of phase: %s' % repr((fq, bc))
            # place the barcode info into the name -- before any whitespace -- delimited by a |
            name, sfx = fq.name.split(' ',1)
            name += '|' + bc[1] + '|' + bc[2]
            if bc[1] == '*':
                handle = reject_handle
            else:
                code = bc[1].split(':')
                lookup = code[1] + code[2] + code[3]
                if lookup in group_map:
                    handle = group_map[lookup]
                elif code[3] != '*':
                    handle = sample_only_handle
                else:
                    handle = reject_handle
            if len(fq.seq) > 0:
                handle.write('@' + name + ' ' + sfx + '\n')
                handle.write(fq.seq + '\n')
                handle.write('+\n')
                handle.write(fq.qual + '\n')
    

if __name__ == '__main__':
    main(get_args())
