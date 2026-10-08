#!/usr/bin/env python3
"""Combine OTU counts, consensus sequences, top-hit, SINTAX, ITSx and chimera
results into the final tables.

Outputs (in --outdir):
  OTU_table.tsv        OTU x sample counts + SINTAX taxonomy (cut-off applied)
  OTU_info.tsv         per OTU: size, number of merged pre-merge OTUs, polishing, length, ITSx regions, chimera
                       flag, top hit (accession, % identity, coverage, taxonomy),
                       SINTAX taxonomy with confidences
  OTU_consensus.fasta  consensus sequences named OTU_n
"""

import argparse
import os
import re
from collections import OrderedDict

RANKS = ['d', 'k', 'p', 'c', 'o', 'f', 'g', 's']
RANK_NAMES = {'d': 'domain', 'k': 'kingdom', 'p': 'phylum', 'c': 'class', 'o': 'order',
              'f': 'family', 'g': 'genus', 's': 'species'}

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--consensus', required=True)
parser.add_argument('--counts', required=True)
parser.add_argument('--top_hit', required=True)
parser.add_argument('--sintax', required=True)
parser.add_argument('--itsx_positions', required=True)
parser.add_argument('--uchime', required=True)
parser.add_argument('--outdir', required=True)
args = parser.parse_args()


def otu_of(label):
    return label.split(';', 1)[0]


def parse_tax(tax_string):
    """'d:Fungi(1.00),p:Ascomycota(0.98)' or 'd:Fungi,p:Ascomycota' -> {rank: (name, conf)}"""
    out = {}
    for part in tax_string.strip().strip(';').split(','):
        m = re.match(r'^([a-z]):(.+?)(?:\(([\d.]+)\))?$', part)
        if m:
            out[m.group(1)] = (m.group(2), m.group(3))
    return out


# consensus sequences + header annotations
cons = OrderedDict()
with open(args.consensus) as fh:
    label = None
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            label = line[1:]
            fields = dict(f.split('=', 1) for f in label.split(';')[1:] if '=' in f)
            cons[otu_of(label)] = {'size': fields.get('size', ''), 'polished': fields.get('polished', ''),
                                   'members': fields.get('members', '1'), 'seq': ''}
        elif label:
            cons[otu_of(label)]['seq'] += line

# counts
with open(args.counts) as fh:
    samples = fh.readline().rstrip('\n').split('\t')[1:]
    counts = {}
    for line in fh:
        f = line.rstrip('\n').split('\t')
        counts[f[0]] = f[1:]

# top hit: query target id alnlen qcov tcov mism opens
top = {}
with open(args.top_hit) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        otu = otu_of(f[0])
        if otu in top:
            continue
        if f[1] == '*':
            top[otu] = {'acc': 'no_hit', 'id': '', 'qcov': '', 'tax': {}}
        else:
            acc, _, tax = f[1].partition(';tax=')
            top[otu] = {'acc': acc, 'id': f[2], 'qcov': f[4], 'tax': parse_tax(tax)}

# SINTAX: query, tax with confidences, strand, tax after cut-off
sintax = {}
with open(args.sintax) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t') + ['', '', '']
        sintax[otu_of(f[0])] = {'full': parse_tax(f[1]), 'cut': parse_tax(f[3])}

# ITSx positions: label, length, SSU: .., ITS1: a-b, 5.8S: .., ITS2: .., LSU: .., [notes]
itsx = {}
with open(args.itsx_positions) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        regions = {}
        for part in f[2:]:
            if ':' in part:
                k, v = part.split(':', 1)
                regions[k.strip()] = v.strip()
        itsx[otu_of(f[0])] = regions

# uchime3: score, query, ..., Y/N/? in last column
chimera = {}
with open(args.uchime) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        chimera[otu_of(f[1])] = f[-1]

used_ranks = [r for r in RANKS if any(r in t['full'] for t in sintax.values()) or any(r in t['tax'] for t in top.values())]

os.makedirs(args.outdir, exist_ok=True)

with open(os.path.join(args.outdir, 'OTU_consensus.fasta'), 'w') as out:
    for otu, c in cons.items():
        out.write(f'>{otu}\n{c["seq"]}\n')

with open(os.path.join(args.outdir, 'OTU_table.tsv'), 'w') as out:
    out.write('\t'.join(['OTU'] + samples + ['total'] + [RANK_NAMES[r] for r in used_ranks]) + '\n')
    for otu in cons:
        row = counts.get(otu, ['0'] * len(samples))
        tax = sintax.get(otu, {}).get('cut', {})
        out.write('\t'.join([otu] + row + [str(sum(map(int, row)))] +
                            [tax.get(r, ('', None))[0] for r in used_ranks]) + '\n')

with open(os.path.join(args.outdir, 'OTU_info.tsv'), 'w') as out:
    header = ['OTU', 'size', 'merged_from', 'polished', 'length', 'ITSx_ITS1', 'ITSx_5.8S', 'ITSx_ITS2', 'ITSx_detected',
              'chimera_uchime3', 'tophit_accession', 'tophit_id', 'tophit_qcov', 'tophit_taxonomy']
    header += [f'sintax_{RANK_NAMES[r]}' for r in used_ranks] + [f'sintax_{RANK_NAMES[r]}_conf' for r in used_ranks]
    out.write('\t'.join(header) + '\n')
    for otu, c in cons.items():
        reg = itsx.get(otu, {})
        detected = 'yes' if any(reg.get(k, 'Not found') != 'Not found' for k in ('ITS1', '5.8S', 'ITS2')) else 'no'
        th = top.get(otu, {'acc': 'no_hit', 'id': '', 'qcov': '', 'tax': {}})
        th_tax = ','.join(f'{r}:{th["tax"][r][0]}' for r in RANKS if r in th['tax'])
        sx = sintax.get(otu, {}).get('full', {})
        row = [otu, c['size'], c['members'], c['polished'], str(len(c['seq'])),
               reg.get('ITS1', ''), reg.get('5.8S', ''), reg.get('ITS2', ''), detected,
               chimera.get(otu, ''), th['acc'], th['id'], th['qcov'], th_tax]
        row += [sx.get(r, ('', None))[0] for r in used_ranks] + [sx.get(r, ('', ''))[1] or '' for r in used_ranks]
        out.write('\t'.join(row) + '\n')

print(f'{len(cons)} OTUs x {len(samples)} samples written to {args.outdir}')
