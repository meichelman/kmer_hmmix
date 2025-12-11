#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1]
runs = []
arc_state = False
run = 0

with open(filename) as f:
    for line in f:
        if line.startswith('start'):
            continue
        start, end, count, mutrate, human_p, arc_p, post_state, viterbi_state = line.strip().split('\t')
        if arc_state:
            if post_state == 'Archaic':
                run += 1
            else:
                runs.append(run)
                arc_state = False
                run = 0
        elif post_state == 'Archaic':
            arc_state = True
            run += 1
        
# If file ends while still in an archaic run, append it
if arc_state and run > 0:
    runs.append(run)

# Plot histogram
plt.hist(runs, bins=range(1, max(runs) + 1), color='skyblue', edgecolor='black', log=True)
plt.xlabel('Length of consecutive Archaic runs (windows)')
plt.ylabel('Frequency')
plt.title('Histogram of Archaic run lengths')
output = filename.replace('probs_and_path.tsv', 'consecutive_arc_states_hist.png')
plt.savefig(output, dpi=300)