# Changelog

## 0.1.1

- Container images are now published to the GitHub Container Registry
  (`ghcr.io/roylejw/poreflow-its-main`, `-medaka`) by a GitHub Actions workflow on every release,
  and the `docker`, `singularity` and `apptainer` profiles use them by default. In 0.1.0 these
  profiles pointed at image names that were never published.
- `ps` added to the images, so Nextflow can record task resource use.
- README: new section on the container images.

## 0.1.0

First version for full-length ITS:

- linked primer trimming with read orientation; optional 3' primer; concatemer removal
- full ITS extracted by cutting at the ITS4 site (replaces barrnap 18S extraction)
- streaming quality sort; quality-ordered clustering (vsearch --cluster_smallmem --usersort)
- batched medaka polishing (single model load, GPU optional)
- merging of near-identical polished consensus OTUs
- ITSx annotation, uchime3 chimera flag, best-hit (usearch_global) and SINTAX taxonomy
- Nextflow implementation: local/SLURM/PBS/SGE/LSF, conda/Docker/Singularity/Apptainer/Wave, scratch execution
