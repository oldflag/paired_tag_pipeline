set -e -x -v
reachtools_refdir=/home/chartl/repos/Paired-Tag/refereces/
reachtools_bin=/home/chartl/repos/Paired-Tag/reachtools/reachtools
genome_ref=/home/share/storages/2T/genome/human/GRCh38.primary_assembly.genome  # has .fa and .bwt extensions
bin_txt=/home/share/storages/2T/genome/human/GRCh38_10kb.txt
transcriptome=/home/share/storages/2T/genome/human/gencode.v39.annotation.txt
transcriptome_index=/home/share/storages/2T/genome/human/star_index

libid=$1
libtype=$2
fq1=$3
fq2=$4

if [[ "${libtype}" != "RNA" ]]; then
  if [[ "${libtype}" != "DNA" ]]; then
    echo "libtype (arg2) must be one of RNA or DNA"
    exit 2
  fi
fi

if [ ! -e "${fq1}" ]; then
  echo "${fq1}" does not exist
  exit 2
fi

if [ ! -e "${fq2}" ]; then
  echo "${fq2}" does not exist
  exit 2
fi

pfx="${libid}_${libtype}"
ln -s $fq1 ${pfx}_R1.fq.gz
ln -s $fq2 ${pfx}_R2.fq.gz

p=$reachtools_refdir $reachtools_bin combine3 ${pfx}
zcat ${pfx}_combined.fq.gz | bowtie -x $reachtools_refdir/cell_id_full_407 - --norc -m 1 -v 1 -S ${pfx}_BC.sam
p=$reachtools_refdir $reachtools_bin convert2 ${pfx}_BC.sam
trim_galore ${pfx}_BC_cov.fq.gz

# map to genome
if [[ "${libtype}" == "DNA" ]]; then
    bowtie2 -x ${genome_ref} -U ${pfx}_BC_cov_trimmed.fq.gz --no-unal -p 8 -S ${pfx}.sam
    samtools sort "${pfx}.sam" -o "${pfx}_sorted.bam" -T ${pfx}_tmpsrt
    $reachtools_bin rmdup2 "${pfx}_sorted.bam"
    $reachtools_bin bam2Mtx2 "${pfx}_sorted_rmdup.bam" "${bin_txt}"
    # cleanup
    rm "${pfx}_combined.fq.gz" "${pfx}_BC.sam" "${pfx}_BC_cov.fq.gz" "${pfx}_BC_cov_trimmed.fq.gz" "${pfx}.sam" "${pfx}_sorted.bam"
elif [[ "${libtype}" == "RNA" ]]; then
     trim_galore -a AAAAAAAAAAAAAAAACCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN ${pfx}_BC_cov_trimmed.fq.gz ### trim oligo-dT primer
     trim_galore -a CCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN ${pfx}_BC_cov_trimmed_trimmed.fq.gz ## trim N6 primer
     STAR  --runThreadN 6 --genomeDir ${transcriptome_index} --readFilesIn ${pfx}_BC_cov_trimmed_trimmed_trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix ${pfx}_ --outSAMtype BAM Unsorted
    samtools view -h -F 256 ${pfx}_Aligned.out.bam -b > ${pfx}_clean.bam
    samtools sort ${pfx}_clean.bam -o ${pfx}_sorted -T ${pfx}_tmpsrt
    $reachtools_bin rmdup2 ${pfx}\_sorted.bam
    $reachtools_bin bam2Mtx2 ${pfx}\_sorted_rmdup.bam ${transcriptome}
    # cleanup
    rm "${pfx}_combined.fq.gz" "${pfx}_BC.sam" "${pfx}_BC_cov.fq.gz" "${pfx}_BC_cov_trimmed.fq.gz"
    rm "${pfx}_BC_cov_trimmed_trimmed.fq.gz" "${pfx}_BC_cov_trimmed_trimmed_trimmed.fq.gz"
    rm "${pfx}_Aligned.out.bam" "${pfx}_clean.bam" "${pfx}_sorted.bam"
else
    echo "Impossible!"
    exit 1
fi
