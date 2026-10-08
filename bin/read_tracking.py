#!/usr/bin/env python3
"""Per-sample read counts through every pipeline step (TSV to stdout).

Usage: read_tracking.py --cluster_summary cluster_summary.tsv --counts <sample>_read_counts.tsv ...
"""

import argparse
import csv

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--cluster_summary', required=True)
parser.add_argument('--counts', nargs='+', required=True, help='per-sample *_read_counts.tsv from preprocessing')
args = parser.parse_args()

rows, steps = {}, []
for path in args.counts:
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            rows.setdefault(r['sample'], {})[r['step']] = r['reads']
            if r['step'] not in steps:
                steps.append(r['step'])

with open(args.cluster_summary) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['sample'] in rows:
            rows[r['sample']]['passed_quality'] = r['reads_clustered']
            rows[r['sample']]['in_OTUs'] = r['reads_in_kept_OTUs']
steps += ['passed_quality', 'in_OTUs']

print('\t'.join(['sample'] + steps + ['pct_raw_in_OTUs']))
for s in sorted(rows):
    raw = int(rows[s].get('raw', 0) or 0)
    kept = int(rows[s].get('in_OTUs', 0) or 0)
    pct = f'{100 * kept / raw:.1f}' if raw else ''
    print('\t'.join([s] + [rows[s].get(st, '0') for st in steps] + [pct]))
