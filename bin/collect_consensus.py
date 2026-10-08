#!/usr/bin/env python3
"""Write one consensus per OTU, in OTU order, taking the medaka-polished sequence
where available and the unpolished centroid otherwise (flagged polished=no).

Usage: collect_consensus.py <medaka_consensus.fasta> <drafts.fasta> <otu_map.tsv> <out.fasta>
"""

import sys


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


polished, drafts = read_fasta(sys.argv[1]), read_fasta(sys.argv[2])
n_fail = 0
with open(sys.argv[3]) as otu_map, open(sys.argv[4], 'w') as out:
    next(otu_map)
    for line in otu_map:
        otu, _cluster, size, _centroid = line.rstrip('\n').split('\t')
        if polished.get(otu):
            seq, flag = polished[otu], 'yes'
        else:
            seq, flag = drafts[otu], 'no'
            n_fail += 1
        out.write(f'>{otu};size={size};polished={flag}\n{seq}\n')
print(f'{len(drafts)} OTUs, {len(drafts) - n_fail} polished, {n_fail} kept unpolished centroid')
