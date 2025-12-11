#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1]
arc_kmers = []
human_kmers = []

with open(filename) as f:
    for line in f:
        if line.startswith('start'):
            continue
        start, end, count, mutrate, human_p, arc_p, post_state, viterbi_state = line.strip().split('\t')
        count = int(count)
        if post_state == 'Archaic':
            arc_kmers.append(count)
        else:
            human_kmers.append(count)

plt.hist(human_kmers,
         bins=range(1, 50),
         color='skyblue',
         edgecolor='black',
         density=True)   # <-- normalize to probability
plt.xlabel('Number of kmers')
plt.ylabel('Probability')
output = filename.replace('probs_and_path.tsv', 'human_kmers_probs.png')
plt.savefig(output, dpi=300)
plt.close()

plt.hist(arc_kmers,
         bins=range(100),
         color='skyblue',
         edgecolor='black',
         density=True)
plt.xlabel('Number of kmers')
plt.ylabel('Probability')
output = filename.replace('probs_and_path.tsv', 'arc_kmers_probs.png')
plt.savefig(output, dpi=300)
plt.close()
