import pandas as pd
from argparse import ArgumentParser
import seaborn as sbn
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = ArgumentParser()
parser.add_argument('rna_metrics', help='The RNASeQC metrics csv file')
parser.add_argument('outpdf', help='The output pdf file')
parser.add_argument('--spikein', help='Distinguish spike-in controls', action='store_true')

args = parser.parse_args()

dat = pd.read_csv(args.rna_metrics)
dat.loc[:, 'library'] = dat.Sample.map(lambda x: x.split('__')[0] if 'unknown' not in x else 'Unassigned')
if args.spikein:
    dat.loc[:, 'library'] = dat.library + '_' + dat.is_spikein
library_order = sorted(dat.library.unique())
dat = dat[['unknown' not in x for x in dat.Sample]]
with PdfPages(args.outpdf) as pdf:
    try:
      fig = sbn.swarmplot(x='library', y='Mapping Rate', hue='library', order=library_order, data=dat)
      fig.get_legend().remove(); plt.xticks(rotation=90);plt.tight_layout();pdf.savefig(); plt.figure()
    except:
      pass
    try:
      fig = sbn.swarmplot(x='library', y='High Quality Intragenic Rate', hue='library', order=library_order, data=dat)
      fig.get_legend().remove(); plt.xticks(rotation=90);plt.tight_layout();pdf.savefig(); plt.figure()
    except:
      pass
    try:
      fig = sbn.swarmplot(x='library', y='High Quality Exonic Rate', hue='library', order=library_order, data=dat)
      fig.get_legend().remove(); plt.xticks(rotation=90);plt.tight_layout();pdf.savefig(); plt.figure()
    except:
      pass
    try:
      fig = sbn.swarmplot(x='library', y='High Quality Intergenic Rate', hue='library', order=library_order, data=dat)
      fig.get_legend().remove(); plt.xticks(rotation=90);plt.tight_layout();pdf.savefig(); plt.figure()
    except:
      pass
    try:
      fig = sbn.swarmplot(x='library', y='Expression Profiling Efficiency', hue='library', order=library_order, data=dat)
      fig.get_legend().remove(); plt.xticks(rotation=90);plt.tight_layout();pdf.savefig(); plt.figure()
    except:
      pass
    try:
      fig = sbn.swarmplot(x='library', y="Median 3' bias", hue='library', order=library_order, data=dat)
      fig.get_legend().remove(); plt.xticks(rotation=90);plt.tight_layout();pdf.savefig()
    except:
      pass
