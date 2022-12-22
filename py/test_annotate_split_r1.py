# Unit tests for annotate_split_r1
import gzip
import numpy as np
import os

import annotate_split_R1 as asr

def create_barcode_file(cps_list, reads_per_cell=10):
    # given a list of cells per sample, create a barcode .csv.gz
    # this file is of the form <readname>,umi:w1:w2:sn,umi:w1b:w2b:sb
    reads = list()
    for i_sample, n_cell in enumerate(cps_list):
        w1, w2 = 1, 1
        for i_cell in range(n_cell):
            w2 += 1
            if w2 > 128:
                w1, w2 = 1 + w1, 1
            for _ in range(reads_per_cell):
                reads.append(f'read_{i_sample}_{len(reads)},AACCTTGG:BC{w1}:BC{w2}:SM{1+i_sample},AACCTTGG:NNNNNNN:NNNNNN:NNNNNN')
    gen = np.random.RandomState(reads_per_cell)
    gen.shuffle(reads)
    tfile = 'tmp_barcode_test.csv.gz'
    out = gzip.open(tfile, 'wt')
    for read in reads:
        out.write(read + '\n')
    return tfile


def test_compute_cells_per_sample():
    """
    Test that `compute_cells_per_sample` constructs the correct
    table of cell counts

    """
    TEST_CELLS_PER_SAMPLE = [
      [1], [125, 1], [12, 14, 18, 2],
      [21, 13, 13, 200, 80],
      [256]*6,
      [256]*6 + [1]*6,
      [1] * 10,
      [1000, 2000]
    ]


    for cell_list in TEST_CELLS_PER_SAMPLE:
        exp_nsample = len(cell_list)
        exp_totcell = sum(cell_list)
        exp_cps = {f'SM{1+x}': y for x, y in enumerate(cell_list)}
        for rps in [1, 10, 100, 500]:
            bcfile = create_barcode_file(cell_list)
            tot_cell, n_sam, d_cps = asr.compute_cells_per_sample(bcfile)
            assert exp_totcell == tot_cell
            assert exp_nsample == n_sam
            for d in exp_cps:
                assert exp_cps[d] == d_cps[d]
            os.system('rm ' + bcfile)


def int2dna(x, l=8):
    s = '' if x > 0 else 'A'
    while x > 0:
        l = x % 4
        s = 'ACGT'[l] + s
        x = int((x-l)/4)
    if len(s) < l:
        s = 'T' * (l - len(s)) + s
    return s
    


def write_barcodes(num_bcs):
    of = 'test_bcs.fa'
    with open(of, 'wt') as out:
        for i in range(num_bcs):
            u = int2dna(i)
            out.write(f'>BC{1+i}\n{u}\n')
    return of


def test_get_group_map():
    """
    Test that `get_group_map` properly splits samples into output fastq files

    """
    SAMPLE_GROUPS = [
       {'SM1': 1,
        'SM2': 1},
       {'SM1': 1,
        'SM2': 2,
        'SM3': 1},
       {'SM1': 4,
        'SM2': 4,
        'SM3': 2},
       {'SM1': 8,
        'SM2': 2,
        'SM3': 4},
       {'SM1': 3,
        'SM2': 7,
        'SM3': 5}]
    N_P_BC = [16, 32, 64, 25, 33, 47, 41, 101]
    with open('digest.csv', 'wt') as out:
        out.write("""assay_id,antibody_name,sequence_library_id,barcode,well
SM0,test,test,GAGG,0
SM1,test,test,AACC,1
SM2,test,test,ATGC,2
SM3,test,test,CATA,3
SM4,test,test,TACC,4
""")
    
    for grouping in SAMPLE_GROUPS:
        for n_barcode in N_P_BC:
            bcf = write_barcodes(n_barcode)
            res, _, _ = asr.get_group_map(bcf, 'digest.csv', grouping, 'test', 'test', open_for_writing=False)
            for sid, exp_nh in grouping.items():
                res_subset = {k: v for k, v in res.items() if k[-3:] == sid}
                revmap = dict()
                for k, v in res_subset.items():
                    if v not in revmap:
                        revmap[v] = list()
                    revmap[v].append(k)
                assert len(revmap) == exp_nh, 'cell mapping size'
                for v in revmap:
                    print((v, len(revmap[v])))
                    if n_barcode % exp_nh == 0:
                        assert len(revmap[v]) == int(n_barcode**2/exp_nh)  # should be an exact match
                    else:
                        l1 = len(revmap[v])
                        l2 = int(n_barcode**2/exp_nh)
                        assert (l1 - l2)/l2 < 0.05  # no more than 5% different
    os.system(f'rm {bcf} digest.csv')


def test_get_base():
    """
    test potential edge cases of the base functionality
    """
    BASES = ['normal', 'thing_with_unders', 'thing withspace',
             'thing_with.dot', 'thing_with_R1bad_info', 
             'double_suffix_R1.fastq.gz']
    SUFFIX = ['_1.fq.gz', '_1.FQ.gz', '_1.fq', '_1.FQ',
              '_1.fastq.gz', '_1.fastq', '_R1.fastq',
              '_R1.fastq.gz', '_R1.fq', '_R1.fq.gz',
              '_R1.FQ.gz', '_R1.fq']
    for base in BASES:
        for sf in SUFFIX:
            newbase = asr.get_base(base + sf)
            assert base == newbase
