"""
This CLT annotates a Paired-Tag R1 fastq.gz with extracted barcodes from R2, and additionally
splits the file into mulitple output fastQs for use with parallel execution of the pipeline.

The output fastQs will be split per-sample and, if there are sufficiently many estimated cells,
then again within sample. The sample ID and antibody will be encoded into the fastQ names.

"""
from argparse import ArgumentParser
from utils import read_fasta, read_fastq, xopen, dxopen
from itertools import product, islice
from csv import DictReader
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
    parser.add_argument('sample_manifest', help='The sample manifest file, containing sample id, barcode, and antibody')
    parser.add_argument('--cells_per_fastq', help='The target # cells / fastQ for output', default=10000, type=int)
    parser.add_argument('--noestimate', help='Do not estimate total # of cells; compute exact numbers', action='store_true')
    parser.add_argument('--sequence_id', help='The sequence id. Defaults to the {this_thing}_R1.fastq.gz', default=None)
    parser.add_argument('--outdir', help='The output directory', default=None)
    parser.add_argument('--min_rlen', help='Minimum R1 read length for output', default=30, type=int)

    return parser.parse_args()



def compute_cells_per_sample(bcfile):
    per_sample_umi = dict()
    with xopen(bcfile, 'rt') as handle:
        for entry in handle:
            bcinfo = entry.split(',')[1]
            if '*' not in bcinfo:
                codes = bcinfo.split(':')                                                                                                                                                                                                                       # umi:bc1:bc2:sm:re
                sample, cell, umi = codes[3], codes[1] + codes[2], codes[0]
                if sample not in per_sample_umi:
                    per_sample_umi[sample] = dict()
                per_sample_umi[sample][cell] = setadd(per_sample_umi[sample].get(cell, set()), codes[0])
    n_samples = len(per_sample_umi)
    cps = {sm: len(per_sample_umi[sm])
           for sm in per_sample_umi}
    total_cells = sum(cps.values())
    return total_cells, n_samples, cps 


def get_group_map(pool_fasta, sample_digest_file, groups_per_sample, out_base, library_id, open_for_writing=True):
    """
    Produce a dictionary of output handles for groups of cell ids

    Input
    ---------
    :pool_fasta: Fasta file of pool barcodes
    :sample_digest_file: csv-formatted digest file containing the fields
                         `assay_id`, `antibody_name`, `sequence_library_id`, `barcode`, `well`
    :groups_per_sample: mapping {sample_id -> n_grp} number of fastq files per
                        sample to write
    :out_base: base for the output file-names
    :library_id: the library id for this run (for matching in the sample digest)
    :open_for_writing: Whether to return open file-handles or just the files. Used for testing

    Output
    ---------
    Map of {bc1:bc2:sm -> handle}, handle to (nomatch) reads, handle to unmatched samples

    """
    bc1_f, bc2_f = list(read_fasta(pool_fasta)), list(read_fasta(pool_fasta))
    nw = len(bc1_f) * len(bc2_f)
    digest_recs = [rec for rec in DictReader(open(sample_digest_file)) 
                   if rec['sequence_library_id'] == library_id]
    out_files, handles = list(), dict()
    for digest_record in digest_recs:
        assay_id, antibody = digest_record['assay_id'], digest_record['antibody_name']
        sample_out = list()
        cpg = int(nw / groups_per_sample.get(digest_record['assay_id'], 1))
        rem = nw % cpg
        if rem > 0 and rem < cpg/2:
            cpg = int(cpg + max(1, rem/groups_per_sample.get(digest_record['assay_id'], 1)))
        rem = nw % cpg
        for i, group_list in enumerate(grouped_product([digest_record['assay_id']], 
                                       bc1_f, bc2_f, 
                                       n=cpg)):
            of = f'{out_base}__{assay_id}__{antibody}__{i+1}.fq.gz'
            hdl = dxopen(of, 'wt') if open_for_writing else of
            out_files.append(of)
            for j, (aid, bc1, bc2) in enumerate(group_list):
                handles[bc1.name + bc2.name + aid] = hdl

    sys.stderr.write(f'Writing to {len(out_files) + 2} fastQ files\n') 
    unknown_barcode = f'{out_base}__unknown__unlinked__1.fq.gz'
    unknown_sample = 'f{out_base}__sample__unlinked__1.fq.gz'
    ubh = dxopen(unknown_barcode, 'wt') if open_for_writing else unknown_barcode
    ush = dxopen(unknown_sample, 'wt') if open_for_writing else unknown_sample

    return handles, ubh, ush
   

def get_base(fn, pfx=None):
    """
    Get the base of a fastq file by stripping:
      .gz .fq .fastq. .FQ .fastQ _R1 _1

    and optionally adding a prefix (dir)
    """
    pfx = '' if pfx is None else pfx.strip('/') + '/'  # ensure just one /
    suffix_list_order = ['.gz', '.fq', '.fastq', '.FQ', 'FASTQ', '_1', '_R1']
    bf = fn
    for sfx in suffix_list_order:
        k = len(sfx)
        if bf[-k:] == sfx:
            bf = bf[:-k]

    return pfx + bf
               

def main(args):
    print(args)
    if args.noestimate:
        ncell, _, cell_per_sam = compute_cells_per_sample(args.bc_csv)
    else:
        ncell, _, cell_per_sam = estimate_via_recapture(args.bc_csv, n=1200000)  # use 1.2M reads to estimate # of cells


    #n_groups = 1 + int(ncell / args.cells_per_fastq)
    groups_per_sample = {sm: 1 + int(nc/args.cells_per_fastq) for sm, nc in cell_per_sam.items()}
    n_groups = sum(groups_per_sample.values())
    sys.stderr.write(f'Using {ncell} estimated total cells and {n_groups} groups\n')

    base = get_base(args.R1_fastq, args.outdir) if args.sequence_id is None else args.outdir + '/' + args.sequence_id
    lib = args.sequence_id if args.sequence_id is not None else base.split('/')[-1]
    print(lib)

    group_map, reject_handle, sample_only_handle = get_group_map(args.well_bc, args.sample_manifest, groups_per_sample, base, lib)

    fq_in = read_fastq(args.R1_fastq)
    with xopen(args.bc_csv, 'rt') as bcinput:
        bc_in = (x.strip().split(',') for x in bcinput)

        for nread, (fq, bc) in enumerate(zip(fq_in, bc_in)):
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
            if len(fq.seq) >= args.min_rlen:
                rname = '@' + name + ' ' + sfx
                if len(rname) > 100:
                    rname = f'read{1+nread}|{bc[1]}|{bc[2]}'
                handle.write('@' + rname + '\n')
                handle.write(fq.seq + '\n')
                handle.write('+\n')
                handle.write(fq.qual + '\n')
    

if __name__ == '__main__':
    main(get_args())
