from contextlib import contextmanager
import gzip
from collections import namedtuple

fasta_record = namedtuple('fasta_record', ('name', 'seq'))
fastq_record = namedtuple('fastq_record', ('name', 'seq', 'qual'))

@contextmanager
def xopen(fn, *args, **kwargs):
    if fn[-2:] == 'gz':
        hdl = gzip.open(fn, *args, **kwargs)
    else:
        hdl = open(fn, *args, **kwargs)

    try:
        yield hdl
    finally:
        hdl.close()


def dxopen(fn, *args, **kwargs):
    if fn[-2:] == 'gz':
        return gzip.open(fn, *args, **kwargs)
    return open(fn, *args, **kwargs)


def read_fasta(fn):
    name, seq = None, ''
    with xopen(fn, 'rt') as f_in:
        for line in f_in:
            if name is None:
                assert line[0] == '>'
                name = line.strip()[1:]
            elif line[0] == '>':
                yield fasta_record(name, seq)
                name, seq = line.strip()[1:], ''
            else:
                seq += line.strip()
    if name is not None:
        yield fasta_record(name, seq)


def read_fastq(fn):
    hdl = read_fastq_inner(fn)
    next(hdl) # ignore the first one
    for rec in hdl:
        yield rec


def read_fastq_inner(fn):
    name, seq, qual = None, '', ''
    with xopen(fn, 'rt') as f_in:
        for i, line in enumerate(f_in):
            if (i % 4) == 0:
                assert line[0] == '@', line
                yield fastq_record(name, seq, qual)
                seq, qual = '', ''
                name = line.strip()[1:]
            elif (i % 4) == 1:
                seq = line.strip()
            elif (i % 4 == 3):
                qual = line.strip()
    yield fastq_record(name, seq, qual)
