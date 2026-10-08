# Changelog

## 0.1.0

First version for full-length ITS:

- linked primer trimming with read orientation; optional 3' primer; concatemer removal
- full ITS extracted by cutting at the ITS4 site (replaces barrnap 18S extraction)
- streaming quality sort; quality-ordered clustering (vsearch --cluster_smallmem --usersort)
- batched medaka polishing (single model load, GPU optional)
- merging of near-identical polished consensus OTUs
- ITSx annotation, uchime3 chimera flag, best-hit (usearch_global) and SINTAX taxonomy
- Nextflow implementation: local/SLURM/PBS/SGE/LSF, conda/Docker/Singularity/Apptainer/Wave, scratch execution
