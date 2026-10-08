#!/bin/bash
# Build a small test dataset for `-profile test` from your own reads and reference:
#   - N random reads from each fastq(.gz) in <reads_dir> (one file per sample)
#   - a reduced reference: references hit by a subset of those reads + 5000 random ones
# Needs seqkit, cutadapt and vsearch (e.g. the envs/main.yml environment).
#
# Usage: bash tests/make_test_data.sh <reads_dir> <reference.fasta> [reads_per_sample=2000]
set -euo pipefail
READS_DIR=$1; DB=$2; N=${3:-2000}
OUT="$(cd "$(dirname "$0")" && pwd)/data"
mkdir -p "$OUT"
echo "sample,reads" > "$OUT/samplesheet.csv"
for f in "$READS_DIR"/*.f*q*; do
    s=$(basename "$f"); s=${s%%.f*q*}
    total=$(seqkit stats -T "$f" | awk -F'\t' 'NR==2{print $4}')
    p=$(awk -v n="$N" -v t="$total" 'BEGIN{p=3*n/t; if(p>1)p=1; print p}')
    seqkit sample -s 11 -p "$p" "$f" -o "$OUT/tmp.fastq.gz"
    seqkit head -n "$N" "$OUT/tmp.fastq.gz" -o "$OUT/$s.fastq.gz"; rm "$OUT/tmp.fastq.gz"
    echo "$s,$s.fastq.gz" >> "$OUT/samplesheet.csv"
    echo "test reads: $OUT/$s.fastq.gz"
done
# references matched by the ITS of 100 reads per sample, plus random background references
# (reads are first cut to the ITS with the default primers / ITS4 site so they align quickly)
for f in "$OUT"/*.fastq.gz; do seqkit head -n 400 "$f" | seqkit fq2fa; done \
    | cutadapt -e 0.15 --rc -g "ACCWGCGGARGGATCATTA;min_overlap=15" --discard-untrimmed - 2>/dev/null \
    | cutadapt -e 0.2 -a "GCATATCAATAAGCGGAGGA;min_overlap=15" --discard-untrimmed -m 300 - 2>/dev/null \
    > "$OUT/query.fa"
vsearch --usearch_global "$OUT/query.fa" --db "$DB" --id 0.8 --strand plus --maxaccepts 5 \
    --maxrejects 32 --threads 8 --userout "$OUT/hits.txt" --userfields target --quiet
cut -f1 "$OUT/hits.txt" | sort -u > "$OUT/ids.txt"
seqkit sample -s 7 -p 0.004 "$DB" | seqkit seq -n -i >> "$OUT/ids.txt"
sort -u "$OUT/ids.txt" -o "$OUT/ids.txt"
seqkit grep -f "$OUT/ids.txt" "$DB" -o "$OUT/test_db.fasta"
rm "$OUT/query.fa" "$OUT/hits.txt" "$OUT/ids.txt"
echo "test reference: $OUT/test_db.fasta ($(grep -c '>' "$OUT/test_db.fasta") sequences)"
echo "Now run: nextflow run main.nf -profile test,conda"
