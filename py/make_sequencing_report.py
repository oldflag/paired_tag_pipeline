import pandas as pd
import os
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# read in the type QC files
files = os.listdir('.')

typeqc_files = [x for x in files if x.endswith('type_qc.csv')]

type_qc = None
for tcf in typeqc_files:
  sequence_id = tcf.split('.')[0]
  lysis_id = sequence_id[-2:]
  expected_type = 'rna' if 'SR' in sequence_id else 'dna'
  dat = pd.read_csv(tcf)
  dat.columns = ['assay_info', 'inferred_type', 'reads']
  dat['sequence_id'] = sequence_id
  dat['lysis_id'] = lysis_id
  dat['expected_type'] = expected_type
  if type_qc is None:
    type_qc  = dat
  else:
    type_qc = pd.concat([type_qc, dat])

alignment_qc = [x for x in files if 'alignment_qc' in x and 'spike' not in x]
alignment_qc = pd.read_csv(alignment_qc[0], header=None)
# SRA2,SRA2__A1679000870943__H3K27me3__ALL_N1_15e6_BM_cells_95p_GFP_Late__1Aligned,rna,A1679000870943,H3K27me3,N barcodes Q20,9624
alignment_qc.columns = ['sequence_id', 'bam_prefix', 'library_type', 'assay_info', 'antibody_name', 'metric', 'value']
alignment_qc['assay_string'] = alignment_qc.bam_prefix.astype(str).map(lambda x: x[:-len('Aligned')] if x.endswith('Aligned') else x)
# hilariously, the alignment qc file is only used to associate the 'assay_info' to the assay string

# read the duplicates files
dup_files = [x for x in files if 'duplication_metrics' in x and 'spikein' not in x]
dup_qc = None
for dfile in dup_files:
    library_type = 'dna' if 'dna' in dfile else 'rna'
    dat = pd.read_csv(dfile)
    # SRA2__A1679000649344__H3K27me3__161Early_10Aligned_tag.bam,61688,19929,32.3061
    dat.columns = ['bam', 'aligned', 'duplicate', 'dup_pct']
    dat['assay_string'] = dat.bam.astype(str).map(lambda x: x.split('Aligned')[0].split('_tag')[0])
    dat['inferred_type'] = library_type
    if dup_qc is None:
      dup_qc = dat
    else:
      dup_qc = pd.concat([dup_qc, dat])


info2assay = dict([(t, t.split('__')[1]) for t in {x for x in dup_qc.assay_string if 'UNK' not in x}])

## preprocess the type QC
type_qc['reads_match'] = type_qc['reads'] * (type_qc['inferred_type'] == type_qc['expected_type']).astype(int)
type_qc['mtype'] = (type_qc['inferred_type'] == type_qc['expected_type']).astype(int)

# Group by sequence_id and assay_info
grouped = type_qc.groupby(['sequence_id', 'assay_info', 'lysis_id'])

# Calculate the required statistics without using agg
typeqc_proc = grouped[['reads', 'reads_match']].sum().reset_index()
typeqc_proc['pct_contam'] = (1 - typeqc_proc['reads_match'] / typeqc_proc['reads']) * 100


dup_qc['assay_info'] = dup_qc.assay_string.map(lambda x: info2assay.get(x, '*'))
dup_qc['sequence_id'] = dup_qc.assay_string.map(lambda x: x.split('__')[0].strip('_UNK'))

combined_qc = typeqc_proc.merge(dup_qc, on=['assay_info', 'sequence_id'], how='inner')
combined_qc['target'] = combined_qc.assay_string.map(lambda x: x.split('__')[2])
combined_qc['sample_name'] = combined_qc.assay_string.map(lambda x: x.split('__')[3])
combined_qc['assay_name'] = combined_qc['sample_name'] + '_' + combined_qc['target']
bad_barcodes = typeqc_proc[typeqc_proc.assay_info == '*']
bad_bc_map = dict(zip(bad_barcodes.sequence_id, bad_barcodes.reads))

valid_qc = combined_qc[combined_qc.assay_info != '*']

## create the dna tables
library_qc = valid_qc.groupby(['sequence_id', 'lysis_id', 'inferred_type'])[['reads', 'reads_match', 'aligned', 'duplicate']].sum().reset_index()
library_qc.columns = ['sequence_id', 'lysis_id', 'type', 'valid_barcodes', 'libtype_match', 'aligned', 'duplicate']
library_qc['invalid_barcodes'] = library_qc.sequence_id.map(bad_bc_map.__getitem__)
library_qc['reads'] = library_qc.invalid_barcodes + library_qc.valid_barcodes
library_qc = library_qc[['sequence_id', 'lysis_id', 'type', 'reads', 'invalid_barcodes', 
                         'valid_barcodes', 'libtype_match', 'aligned', 'duplicate']]
library_qc['% Valid'] = round(1000 * library_qc.valid_barcodes / library_qc.reads)/10
library_qc['% Contam'] = round( 1000 * (1 - library_qc.libtype_match/library_qc.valid_barcodes))/10
library_qc['% Aligned'] = round(1000 * library_qc.aligned/library_qc.libtype_match)/10
library_qc['% Duplicate'] = round(1000 * library_qc.duplicate/library_qc.aligned)/10
library_qc.sort_values(by='sequence_id').to_csv('LIBRARY_REPORT.csv', index=False)

## create the sample tables
for ainf in combined_qc.assay_info.unique():
    qc_info = combined_qc[combined_qc.assay_info == ainf]
    samid = qc_info.assay_name.values[0]
    sample_qc = qc_info[['sequence_id', 'lysis_id', 'inferred_type', 'reads', 'reads_match', 'aligned', 'duplicate']].copy()
    sample_qc.columns = ['sequence_id', 'lysis_id', 'type', 'valid_barcodes', 'libtype_match', 'aligned', 'duplicate']
    sample_qc['% Contam'] = round(1000 * (1 - sample_qc.libtype_match / sample_qc.valid_barcodes))/10
    sample_qc['% Aligned'] = round(1000 * sample_qc.aligned / sample_qc.libtype_match)/10
    sample_qc['% Duplicate'] = round(1000 * sample_qc.duplicate / sample_qc.aligned)/10
    sample_qc.sort_values(by='sequence_id').to_csv(f'SAMPLE_REPORT.{samid}.csv', index=False)


# Set up LaTeX text rendering
#plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Create a figure and axis
fig, ax = plt.subplots(figsize=(12, 6))

df = library_qc[['sequence_id', 'lysis_id', 'type', 'reads', '% Valid', '% Contam', '% Aligned', '% Duplicate']]
df.columns = ['Library', 'Lysis', 'Type', '# Reads', '% Valid', '% Contam.', '% Aligned', '% Duplicate']
# Plot the table
ccolors = ([['#f0f0f0'] * len(df.columns)] + [['#d3d3d3'] * len(df.columns)]) * int(len(df)/2)  # Alternate-row coloring
print(len(ccolors))
table = ax.table(cellText=df.values,
                 colLabels=df.columns,
                 cellLoc='center',
                 loc='center',
                 cellColours=ccolors,
                 colColours=['#f9ec9a'] * len(df.columns))  # Header row color

# Format table cells
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1.2, 1.2)  # Adjust the scaling as needed for proper cell size

# Remove axis
ax.axis('off')

# Save the figure as a PDF
with PdfPages('LIBRARY_REPORT.pdf') as pdf:
    pdf.savefig(fig, bbox_inches='tight')
