#!/usr/bin/env python

import json
from itertools import groupby


json_file = "$json"
alignments = "$msa"
# standardized workflow file name
output_prefix = alignments.rpartition(".msa.fasta")[0]
regions = set()


def read_fasta(fh):
    for header, group in groupby(fh, lambda line: line[0] == ">"):
        if header:
            line = next(group)
            name = line.strip()
        else:
            seq = "".join(line.strip() for line in group)
            yield name, seq


with open(json_file) as fh:
    gard = json.load(fh)
    total_len = int(gard["input"]["number of sites"])

    for i in range(len(gard["improvements"])):
        start = 0
        for bp in gard["improvements"][f"{i}"]["breakpoints"]:
            end = bp[0]
            regions.add(f"{start}_{end}")
            start = int(end) + 1
        regions.add(f"{start}_{total_len}")

# grab the sequences
seqs = dict()
with open(alignments) as fh:
    for name, seq in read_fasta(fh):
        seqs[name] = seq.replace("u", "t").replace("U", "T")

for region in regions:
    filename = f"{output_prefix}_{region}.fa"
    start, end = region.split("_")
    start = int(start)
    end = int(end)
    with open(filename, "w") as fh:
        for name, seq in seqs.items():
            print(name, seq[start:end+1], sep="\\n", file=fh)
