#!/usr/bin/env python3
"""Turn a vsearch .uc clustering into per-OTU read files for consensus polishing,
plus an OTU x sample count table.

Read labels must carry ';sample=<name>'.

For every cluster with >= --min_size reads (numbered OTU_1.. by decreasing size):
  <out>/clusters/OTU_n.fastq            up to --max_reads member reads (centroid
                                        first, then in input = quality order)
  <out>/clusters/OTU_n_centroid.fasta   the centroid, renamed OTU_n
Tables:
  <out>/otu_map.tsv      OTU, vsearch cluster number, size, centroid label
  <out>/otu_counts.tsv   OTU x sample read counts (all members, not subsampled)
  <out>/cluster_summary.tsv  reads in kept vs. discarded clusters, per sample
"""

import argparse
import os
from collections import defaultdict

from Bio import SeqIO

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-c', '--uc', required=True)
parser.add_argument('-f', '--fastq', required=True, help='uncompressed FASTQ that was clustered')
parser.add_argument('-o', '--outdir', required=True)
parser.add_argument('-s', '--min_size', type=int, default=5)
parser.add_argument('-m', '--max_reads', type=int, default=1000)
args = parser.parse_args()


def read_id(label):
    return label.split(';', 1)[0]


def sample_of(label):
    for field in label.split(';'):
        if field.startswith('sample='):
            return field[len('sample='):]
    raise ValueError(f'no ;sample= annotation in read label: {label}')


centroid = {}                 # cluster number -> centroid label
members = defaultdict(list)   # cluster number -> labels (centroid first)
with open(args.uc) as uc:
    for line in uc:
        f = line.rstrip('\n').split('\t')
        if f[0] == 'S':
            centroid[f[1]] = f[8]
            members[f[1]].insert(0, f[8])
        elif f[0] == 'H':
            members[f[1]].append(f[8])

clusters = sorted(members, key=lambda c: (-len(members[c]), int(c)))
kept = [c for c in clusters if len(members[c]) >= args.min_size]
otu_name = {c: f'OTU_{i}' for i, c in enumerate(kept, 1)}

samples = sorted({sample_of(l) for c in clusters for l in members[c]})
reads_total = defaultdict(int)
reads_kept = defaultdict(int)
counts = {c: defaultdict(int) for c in kept}
for c in clusters:
    for label in members[c]:
        s = sample_of(label)
        reads_total[s] += 1
        if c in counts:
            counts[c][s] += 1
            reads_kept[s] += 1

os.makedirs(os.path.join(args.outdir, 'clusters'), exist_ok=True)

with open(os.path.join(args.outdir, 'otu_map.tsv'), 'w') as out:
    out.write('OTU\tvsearch_cluster\tsize\tcentroid\n')
    for c in kept:
        out.write(f'{otu_name[c]}\t{c}\t{len(members[c])}\t{centroid[c]}\n')

with open(os.path.join(args.outdir, 'otu_counts.tsv'), 'w') as out:
    out.write('OTU\t' + '\t'.join(samples) + '\n')
    for c in kept:
        out.write(otu_name[c] + '\t' + '\t'.join(str(counts[c][s]) for s in samples) + '\n')

with open(os.path.join(args.outdir, 'cluster_summary.tsv'), 'w') as out:
    out.write('sample\treads_clustered\treads_in_kept_OTUs\treads_in_small_clusters\n')
    for s in samples:
        out.write(f'{s}\t{reads_total[s]}\t{reads_kept[s]}\t{reads_total[s] - reads_kept[s]}\n')

index = SeqIO.index(args.fastq, 'fastq', key_function=read_id)
for c in kept:
    name = otu_name[c]
    with open(os.path.join(args.outdir, 'clusters', name + '.fastq'), 'w') as out:
        for label in members[c][:args.max_reads]:
            rec = index[read_id(label)]
            rec.id, rec.description = read_id(label), ''
            out.write(rec.format('fastq'))
    with open(os.path.join(args.outdir, 'clusters', name + '_centroid.fasta'), 'w') as out:
        rec = index[read_id(centroid[c])]
        out.write(f'>{name}\n{rec.seq}\n')
index.close()

print(f'{len(clusters)} clusters, {len(kept)} with >= {args.min_size} reads written to {args.outdir}/clusters')
