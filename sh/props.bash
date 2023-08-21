fastq=$1

tcount=0
declare -A counts
zcat "${fastq}" | grep ^@ | tr '|' '\t' | cut -f2 | tr ':' '\t' | cut -f 5 | sort | uniq -c > tmp
while read count type; do
  counts["${type}"]=$count
  tcount=$(($tcount + $count))
done < tmp

out=""
for key in "${!counts[@]}"; do
  val=$(echo "scale=2; 100*${counts[$key]}/$tcount" | bc)
  out+=" ${key}=${val}%"
done

echo "${fastq}:${out}"
