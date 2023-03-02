import sys
import numpy
import scipy
import scipy.sparse
import scipy.io
import gzip
import numpy as np

# convert a long form: <cell> <gene_list> <count> into cell_id.txt gene_id.txt counts.mtx format

in_txt = sys.argv[1]
if len(sys.argv) > 2:
    feature_type = sys.argv[2]
else:
    feature_type = 'RNA_expression'

genes, cells = list(), list()
geneset, cellset = dict(), dict()
data, i_, j_ = list(), list(), list()

if in_txt[-2:] == 'gz':
    inf = gzip.open(in_txt, 'rt')
else:
    inf = open(in_txt, 'rt')

next(inf)  # header
for n, line in enumerate(inf):
    fields = line.strip().split('\t')
    features = fields[0].split(',')
    count = int(fields[-1])
    dc = count#/len(features)
    if fields[1] not in cellset:
        cells.append(fields[1])
        cellset[fields[1]] = len(cells) - 1
        ic = len(cells)-1
    else:
        ic = cellset[fields[1]]
    for f in features:
        j_.append(ic)
        if f not in geneset:
            genes.append(f)
            geneset[f] = len(genes)-1
            i_.append(len(genes)-1)
        else:
            i_.append(geneset[f])
        data.append(dc)
    if (n+1) % 1000000 == 0:
        print('... ' + str(n+1))


counts = scipy.sparse.coo_matrix((data, (i_, j_)), shape=(len(genes), len(cells)), dtype=float)
counts = scipy.sparse.csc_matrix(counts)

print(counts.shape)
print(len(cells))
print(len(genes))

base = sys.argv[1][:-len('.txt.gz')]

with open('barcodes.txt', 'wt') as out:
    out.write('\n'.join(cells))

with open('genes.txt', 'wt') as out:
    for gene in genes:
        out.write('%s\t%s\t%s\n' % (gene, gene, feature_type))

scipy.io.mmwrite('counts.mtx', counts)

