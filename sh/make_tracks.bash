bamlist=$1
bam2fragpy=$2
outwig=$3
fai=$4

unsorted=$(echo "${outwig}" | sed 's/.bw$/.unsorted.bed/g')
sorted=$(echo "${outwig}" | sed 's/.bw$/.sorted.bed/g')
bedgraph=$(echo "${outwig}" | sed 's/.bw$/.bedgraph/g')
while read bamfile; do
    fragfile=$(echo "${bamfile}" | sed 's/.bam$/.tsv.gz/g')
    bedfile=$(echo "${bamfile}" | sed 's/.bam$/.bed/g')
    python "${bam2fragpy}" "${bamfile}" "${fragfile}" --min_count 20 --nocb
    zcat "${fragfile}" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t.\t."}' | grep "^chr" >> "${unsorted}"
done < "${bamlist}"

cat "${fai}" | cut -f1-2 > genome

bedtools sort -i "${unsorted}" > "${sorted}"
bedtools genomecov -i "${sorted}" -g genome -bg > "${bedgraph}"
bedGraphToBigWig "${bedgraph}" genome "${outwig}"

