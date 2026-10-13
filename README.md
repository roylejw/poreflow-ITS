# poreflow-ITS

**Full-length ITS (ITS1-5.8S-ITS2) metabarcoding of Oxford Nanopore amplicon reads.**

`poreflow-ITS` turns raw or demultiplexed ONT reads from long ITS amplicons into a polished,
taxonomically annotated OTU table. It is a Nextflow pipeline that runs on a laptop or on SLURM, PBS Pro,
SGE or LSF clusters, with conda, Docker, Singularity/Apptainer or automatically built containers.

```
raw reads (per sample: one fastq(.gz) or a folder of dorado/MinKNOW split files)
 │ 1  PREPROCESS     cutadapt: linked primers F...rc(R) with --rc → trimmed, all reads oriented F→R
 │                   discard reads with a second full primer (concatemers)
 │                   cut at the ITS4 site (~40 bp into LSU) → full ITS only
 ▼
 │ 2  CLUSTER        pool samples; sort reads by mean quality; drop < Q20
 │                   vsearch --cluster_smallmem --usersort (best reads become centroids)
 │                   OTUs = clusters ≥ 5 reads
 ▼
 │ 3  POLISH         medaka consensus of every OTU (one batched GPU/CPU run)
 ▼
 │ 4  MERGE_OTUS     merge near-identical polished consensus OTUs (vsearch --cluster_size)
 │    ITSX           ITS1 / 5.8S / ITS2 boundaries
 │    CHIMERA        vsearch --uchime3_denovo
 │    TOPHIT         vsearch --usearch_global best hit (≥ 90 % query coverage),
 │                   in parallel chunks of tophit_chunk_size OTUs
 │    SINTAX         vsearch --sintax with bootstrap confidence
 ▼
 OTU_TABLE           OTU × sample counts, taxonomy, QC flags, read tracking
```

## Quick start

1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) (≥ 24.10, Java 17+) and one of
   conda/mamba, Docker, Singularity or Apptainer.
2. Get a SINTAX-formatted ITS reference, for example
   [EUKARYOME](https://eukaryome.org) `SINTAX_EUK_ITS_v*.fasta` (all eukaryotes) or a
   [UNITE](https://unite.ut.ee) SINTAX release (fungi).
3. Describe your samples in a samplesheet (or just point `--input` at a folder of reads, see [Input](#input)):

   ```csv
   sample,reads
   site1,/data/run1/fastq_pass/barcode01
   site2,/data/run1/fastq_pass/barcode02
   pos_control,/data/run1/pos_control.fastq.gz
   ```

4. Copy `params.example.yaml` to `my_params.yaml` and set `input`, `db` and `outdir`, plus your primers
   if they differ from the defaults.
5. Run with one **executor** profile and one **software** profile:

   ```bash
   # laptop / workstation
   nextflow run main.nf -profile local,conda -params-file my_params.yaml

   # SLURM cluster, GPU polishing, site settings in my_site.config
   nextflow run main.nf -profile slurm,conda -params-file my_params.yaml -c my_site.config
   ```

   Add `-resume` to continue a run after an interruption or after changing late-stage settings.

## Profiles

| Executor | | Software | |
|---|---|---|---|
| `local` | this machine | `conda` / `mamba` | builds `envs/main.yml` and `envs/medaka.yml` on first use |
| `slurm` | SLURM | `docker` | images pulled from `ghcr.io/roylejw/poreflow-its-*` |
| `pbs` | PBS Pro | `singularity` / `apptainer` | the same images, converted automatically on first use |
| `sge` | SGE / UGE (`smp` PE) | `wave` | containers built on the fly from `envs/*.yml` by Seqera Wave (needs internet) |
| `lsf` | LSF | | |

On every cluster profile each task runs in node-local scratch (`$TMPDIR`, `process.scratch = true`) and only
its outputs are copied back to the work directory. The work directory itself (`-work-dir`, or `workDir` in a site
config) must be on storage shared by all compute nodes.

### Containers

The two software images (`main` and `medaka`) are built from `envs/*.yml` by `containers/Dockerfile.*` and published
to the GitHub Container Registry. Docker, Singularity and Apptainer pull them automatically; Singularity/Apptainer convert them to `.sif` images on first
use (set `NXF_SINGULARITY_CACHEDIR` to a shared folder so this happens once). To build them yourself instead, run
`bash containers/build.sh` and point `--container_main` / `--container_medaka` at your tags.

### Site configuration

Cluster-specific settings go in a config file passed with `-c`. Use `conf/site_example.config` as a
template: maximum node size, queue/account, GPU request, job submission rate, shared conda cache, etc.

```groovy
params {
    cluster_options     = '--account=myproject'  // every job
    gpu                 = true
    gpu_cluster_options = '--gpus-per-node=1'    // GPU job only; or '--gres=gpu:1', or PBS '-l select=1:ngpus=1'
}
process {
    resourceLimits = [cpus: 32, memory: 256.GB, time: 72.h]   // largest job allowed
}
```

(`--max_cpus`, `--max_memory` and `--max_time` set the same caps from the command line or a params file.)

### GPU

Only the polishing step (medaka) benefits from a GPU. With `--gpu true` it receives
`gpu_cluster_options` on HPC executors, `--gpus all` with Docker, and `--nv` with Singularity/Apptainer.
For the conda profile, create the medaka environment where CUDA is visible, or set
`CONDA_OVERRIDE_CUDA=12.4` when it is built on a login node without GPUs. The polishing log
(`03_polish/medaka_device.log`) records whether a GPU was used.

## Input

`--input` accepts:

- **a samplesheet** (CSV or TSV) with columns `sample,reads`. `reads` is a fastq(.gz) file or a folder of
  fastq(.gz) files, for example the split files of one barcode from MinKNOW/dorado (`fastq_pass/barcode01`).
  Relative paths are relative to the samplesheet.
- **a folder** containing one fastq(.gz) per sample (sample = file name) and/or one sub-folder of fastq files
  per sample (sample = folder name).

Reads can be raw (adapters and barcodes still attached) or already demultiplexed/trimmed. Everything outside the
primers is removed. Sample names may contain letters, digits, `.`, `_` and `-`.

## Parameters

All parameters can be set in the params file (`-params-file`) or on the command line (`--name value`).

| Parameter | Default | Description |
|---|---|---|
| `input` | – | samplesheet or folder of reads |
| `db` | – | SINTAX-formatted reference FASTA (`>id;tax=d:...,p:...,g:...,s:...;`) |
| `outdir` | `results` | output folder |
| `primer_fwd` | `ACCWGCGGARGGATCATTA` | forward primer (SSU end) |
| `primer_rev` | `TCCTGAGGGAAACTTCG` | reverse primer (LSU) |
| `its_end_site` | `GCATATCAATAAGCGGAGGA` | 3′ cut site defining the end of the ITS (reverse complement of ITS4) |
| `primer_error` | 0.15 | maximum primer error rate |
| `require_rev_primer` | false | also require the 3′ primer (otherwise trimmed only if present) |
| `min_amplicon_len` / `max_amplicon_len` | 1000 / 2500 | length window after primer trimming |
| `min_its_len` / `max_its_len` | 300 / 1500 | length window of the extracted full ITS |
| `min_read_q` | 10 | minimum mean read quality |
| `cluster_id` | 0.97 | read clustering identity |
| `min_cluster_size` | 5 | minimum reads per OTU |
| `max_reads_per_cluster` | 1000 | highest-quality reads per OTU used for polishing |
| `medaka_model` | `r1041_e82_400bps_hac_v5.2.0` | medaka model matching the basecaller (`medaka tools list_models`) |
| `merge_id` | 0.97 | identity for merging polished consensus OTUs (1 = no merging) |
| `tax_query_cov` | 0.9 | minimum query coverage of the top hit |
| `tophit_chunk_size` | 1000 | OTUs per best-hit search job; jobs run in parallel and their hits are joined |
| `sintax_cutoff` | 0.8 | SINTAX bootstrap cut-off for the OTU table taxonomy |
| `sintax_seed` | 1 | random seed for SINTAX, so repeated runs give identical confidences |
| `itsx_taxa` | `all` | ITSx profile set (`-t`), e.g. `F` for fungi only |
| `per_read_tax` | false | also classify every read with SINTAX |
| `cluster_options` | `''` | scheduler options added to every job (e.g. account) |
| `gpu` / `gpu_cluster_options` | false / `''` | GPU polishing and the scheduler flags that request a GPU |
| `max_cpus` / `max_memory` / `max_time` | 16 / 64.GB / 48.h | caps on any single task |
| `conda_main` / `conda_medaka` | – | use existing conda environments instead of building them |
| `container_main` / `container_medaka` | `ghcr.io/roylejw/poreflow-its-*:0.2.0` | container images |

### Adapting to other primers

The defaults target an amplicon from the 3′ end of the SSU to the LSU (~1.7–1.9 kb), from which the
full ITS is cut at the ITS4 site. For another primer pair, set `primer_fwd`, `primer_rev` and
the length windows. If the reverse primer lies upstream of the ITS4 site (e.g. ITS4 itself), set
`its_end_site` to a conserved site just inside the amplicon, or to the reverse complement of the
reverse primer with `min_amplicon_len` lowered accordingly.

## Output

```
results/
├── results/
│   ├── OTU_table.tsv            OTU × sample counts + SINTAX taxonomy (≥ sintax_cutoff)
│   ├── OTU_info.tsv             per OTU: size, merged OTUs, polished?, length, ITSx regions,
│   │                            chimera flag, top hit (accession, % id, coverage, taxonomy),
│   │                            SINTAX names and confidences per rank
│   ├── OTU_consensus.fasta      polished consensus sequences
│   ├── OTU_merge_membership.tsv read-level OTUs merged into each final OTU
│   └── read_tracking.tsv        reads per sample after every step
├── 01_preprocess/<sample>/      full-ITS reads, cutadapt logs, read counts
├── 02_cluster/                  clustering (.uc), OTU map and counts, read qualities
├── 03_polish/                   consensus sequences, medaka logs
├── 04_taxonomy/                 merging, ITSx, uchime3, top hit, SINTAX outputs
└── pipeline_info/               software_versions.yml, execution report, timeline, trace
```

Suggested to keep OTUs with `ITSx_detected == yes` and `chimera_uchime3 != Y`.

## Testing

```bash
bash tests/make_test_data.sh <folder_with_reads> <reference.fasta>   # small subsample + reduced reference
nextflow run main.nf -profile test,conda
```

## Method notes

- **Primer trimming and orientation.** A linked adapter `F...rc(R)` with `--rc` removes everything
  outside the primers and orients every read F→R. The 5′ primer is required. The 3′ primer is optional
  by default, because pre-trimmed or truncated reads often keep only part of it. Reads that still
  contain a second full primer (concatemers) are discarded.
- **Full ITS extraction.** Per-read ITSx is too slow for millions of nanopore reads, so reads are cut at the
  conserved ITS4 site. This removes ~1 kb of conserved LSU that would otherwise dominate clustering identity.
  ITSx is run on the polished consensus sequences to annotate and validate the ITS.
- **Quality-ordered clustering.** Reads are sorted by mean quality, so the most accurate reads become
  centroids (`--usersort`).
- **Clustering run time.** All reads are pooled and clustered together, so OTUs are consistent across
  samples. The run time of this step grows faster than the number of reads (about 1 h for 2 million ITS reads
  on 14 CPUs), and it uses a single job. For very large runs (tens of millions of reads), give CLUSTER
  more time and CPUs in your site config, e.g.
  `process { withName: 'CLUSTER' { cpus = 40; time = 168.h } }` with `max_time` raised to match.
- **Consensus merging.** At nanopore error rates, reads with more error than 1 − `cluster_id` relative to
  the centroid start new clusters. One template is therefore split into several OTUs whose polished
  consensus sequences differ by a few indels. These are merged at `merge_id`, and their counts summed.

## Credits and licence

Distributed under GPL-3.0 (see `LICENSE`).

Please also cite the tools the pipeline uses (versions are recorded in
`pipeline_info/software_versions.yml`): Nextflow, cutadapt, VSEARCH, SeqKit, medaka, minimap2,
SAMtools, ITSx, HMMER, Biopython, and the reference database you used.

Citation for poreflow-ITS: *in preparation* (see `CITATION.cff`).
