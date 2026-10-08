#!/usr/bin/env python3
"""Merge polished OTU consensus sequences that were re-clustered with
vsearch --cluster_size (abundance order), summing their read counts.

Read-level clustering of nanopore reads splits one template into many OTUs
(reads with >1-(id) error vs. the centroid start new clusters); after polishing,
those consensus sequences are near-identical and are merged here. The merged OTU
takes the sequence of its most abundant member.

Outputs (in --outdir):
  consensus.fasta     >OTU_n;size=<total reads>;polished=<yes/no>;members=<k>
  otu_counts.tsv      OTU x sample counts
  otu_membership.tsv  pre-merge OTU -> merged OTU
"""

import argparse
import os
from collections import defaultdict

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--consensus', required=True, help='polished consensus FASTA (>OTU_n;size=..;polished=..)')
parser.add_argument('--counts', required=True, help='pre-merge OTU x sample counts')
parser.add_argument('--uc', required=True, help='vsearch --cluster_size .uc of the consensus sequences')
parser.add_argument('--outdir', required=True)
args = parser.parse_args()


def otu_of(label):
    return label.split(';', 1)[0]


seqs, flags, order = {}, {}, []
with open(args.consensus) as fh:
    name = None
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            name = otu_of(line[1:])
            fields = dict(f.split('=', 1) for f in line[1:].split(';')[1:] if '=' in f)
            flags[name] = fields.get('polished', '')
            seqs[name] = ''
            order.append(name)
        elif name:
            seqs[name] += line

with open(args.counts) as fh:
    samples = fh.readline().rstrip('\n').split('\t')[1:]
    counts = {}
    for line in fh:
        f = line.rstrip('\n').split('\t')
        counts[f[0]] = [int(x) for x in f[1:]]

parent = {}
with open(args.uc) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if f[0] == 'S':
            parent[otu_of(f[8])] = otu_of(f[8])
        elif f[0] == 'H':
            parent[otu_of(f[8])] = otu_of(f[9])

groups = defaultdict(list)
for otu in order:
    groups[parent.get(otu, otu)].append(otu)

totals = {}
for rep, members in groups.items():
    totals[rep] = [sum(counts.get(m, [0] * len(samples))[i] for m in members) for i in range(len(samples))]
reps = sorted(groups, key=lambda r: (-sum(totals[r]), order.index(r)))
new_name = {r: f'OTU_{i}' for i, r in enumerate(reps, 1)}

os.makedirs(args.outdir, exist_ok=True)
with open(os.path.join(args.outdir, 'consensus.fasta'), 'w') as out:
    for r in reps:
        out.write(f'>{new_name[r]};size={sum(totals[r])};polished={flags[r]};members={len(groups[r])}\n{seqs[r]}\n')
with open(os.path.join(args.outdir, 'otu_counts.tsv'), 'w') as out:
    out.write('OTU\t' + '\t'.join(samples) + '\n')
    for r in reps:
        out.write(new_name[r] + '\t' + '\t'.join(map(str, totals[r])) + '\n')
with open(os.path.join(args.outdir, 'otu_membership.tsv'), 'w') as out:
    out.write('pre_merge_OTU\tmerged_OTU\tpre_merge_reads\tis_representative\n')
    for r in reps:
        for m in groups[r]:
            out.write(f'{m}\t{new_name[r]}\t{sum(counts.get(m, []))}\t{"yes" if m == r else "no"}\n')

print(f'{len(order)} consensus OTUs merged into {len(reps)} OTUs')
