# Test data

`make_test_data.sh` builds a small dataset from your own reads, because the
pipeline does not ship sequencing data:

```bash
bash tests/make_test_data.sh <folder_with_one_fastq_per_sample> <SINTAX_reference.fasta> [reads_per_sample]
nextflow run main.nf -profile test,conda
```

It writes `tests/data/` (ignored by git): a random subsample of each file, a samplesheet, and a
reduced reference made of the references hit by a subset of the reads plus ~5000 random ones.
