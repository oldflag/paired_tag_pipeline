#!/bin/bash
if [ -z ${CONDAPATH+x} ]; then
  CONDAPATH="/home/chartl/miniconda3/bin"
  CONDAENV="splitbam"
fi

if [ -z ${REPO+x} ]; then
  REPO="/home/chartl/repos/pipelines/py"
fi

echo source "${CONDAPATH}/activate" "${CONDAENV}"
source "${CONDAPATH}/activate" "${CONDAENV}"
conda activate "${CONDAENV}"

env | grep PATH

python "${REPO}/create_cluster_bam.py" "$1" "$2" "$3"
for bbase in `ls $3/*.bam`; do
  bname="${bbase%%.*}"
  BAMscale scale --bam "${3}/${bbase}" -t 4 --smoothen 20
done
