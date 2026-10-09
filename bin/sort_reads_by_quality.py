#!/usr/bin/env python3
"""Sort a FASTQ by decreasing mean read quality and drop reads below a threshold.

Mean quality is -10*log10(mean per-base error probability). Reads are not all
held in memory: qualities are computed in one pass, then the records are written
in sorted order from a SeqIO index (input must be uncompressed).

Usage: sort_reads_by_quality.py -i in.fastq -o sorted.fastq -q 10 [-t quality.tsv]
"""

import argparse
import math

import numpy as np

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i', '--input_fastq', required=True, help='uncompressed FASTQ')
parser.add_argument('-o', '--output_fastq', required=True)
parser.add_argument('-q', '--min_qual', type=float, default=0)
parser.add_argument('-t', '--quality_table', help='optional TSV: read_id, mean_quality, mean_error_prob')
args = parser.parse_args()

# phred+33 byte value -> error probability
ERR = np.array([10 ** (-(b - 33) / 10) if b >= 33 else 1.0 for b in range(256)])

records = []
with open(args.input_fastq) as handle:
    for title, _seq, qual in FastqGeneralIterator(handle):
        err = float(ERR[np.frombuffer(qual.encode('ascii'), dtype=np.uint8)].mean())
        records.append((-10 * math.log10(err), title.split(None, 1)[0], err))

records.sort(key=lambda x: x[0], reverse=True)

index = SeqIO.index(args.input_fastq, 'fastq')
kept = 0
with open(args.output_fastq, 'wb') as out:
    for q, read_id, _err in records:
        if q < args.min_qual:
            break
        out.write(index.get_raw(read_id))
        kept += 1
index.close()

if args.quality_table:
    with open(args.quality_table, 'w') as out:
        out.write('read_id\tmean_quality\tmean_error_prob\n')
        for q, read_id, err in records:
            out.write(f'{read_id}\t{q:.3f}\t{err:.6g}\n')

print(f'{kept} of {len(records)} reads kept (mean Q >= {args.min_qual}), written to {args.output_fastq}')
