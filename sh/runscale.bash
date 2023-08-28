set -e
set -v

pseudobulk_bam=$1
valid_barcodes=$2

filtered_bam=$(echo "${pseudobulk_bam}" | sed 's/.bam$/.validBC.bam/g')
out_bw=$(echo "${pseudobulk_bam}" | sed 's/.bam$/.validBC.bw/g')

# first : match read names to the valid barcodes
qryf=$(echo "${valid_barcodes}" | sed 's/.txt$/.csv/g' | sed 's/.csv$/.query.csv/g')

echo "Query file: ${qryf}  filtered bam: ${filtered_bam}"

echo @SQ > "${qryf}"
echo @PG >> "${qryf}"
echo @HD >> "${qryf}"

cat "${valid_barcodes}" >> "${qryf}"

samtools view -h "${pseudobulk_bam}" | grep -f "${qryf}" | samtools view -hb > "${filtered_bam}"
samtools index "${filtered_bam}"

BAMscale scale -i "${filtered_bam}" -k custom -F 0.1 -t 6 

