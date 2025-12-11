#!/usr/bin/env python

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from upsetplot import from_memberships, plot


def parse_file(filename):
    arc_prob = []
    kmer_counts = []
    starts = []
    ends = []
    post_states = []
    
    with open(filename) as f:
        
        for line in f:
            if line.startswith('start'):
                continue
            
            if len(line.strip().split('\t')) == 8:
                start, end, count, mutrate, human_p, arc_p, post_state, viterbi_state = line.strip().split('\t')
                arc_prob.append(float(arc_p))
                kmer_counts.append(int(count))
                starts.append(int(start))
                ends.append(int(end))
                post_states.append(post_state)
            elif len(line.strip().split('\t')) <= 5:
                count = line.strip().split('\t')[3]
                kmer_counts.append(int(count))

    return (np.array(arc_prob),
            np.array(kmer_counts),
            np.array(starts),
            np.array(ends),
            np.array(post_states)
           )


def main():
    
    args = sys.argv
    if len(args) < 2:
        sys.exit('\n\nMust input filename\n\n')
    filename = args[1]
    archaics = { 'denisova': [],
                'altai': [],
                'chag': [],
                'vindija': []
               }
    for arc in archaics.keys():
        arc_prob, kmer_counts, starts, ends, post_state = parse_file(filename)
        archaics[arc] = kmer_counts
        # archaics[arc] = post_state
    
    combinations = [[ 'denisova', 'altai', 'chag', 'vindija' ],
         [ 'denisova', 'altai', 'chag' ],
         [ 'denisova', 'altai', 'vindija' ],
         [ 'denisova', 'chag', 'vindija' ],
         [ 'altai', 'chag', 'vindija' ],
         [ 'denisova', 'altai' ],
         [ 'denisova', 'chag' ],
         [ 'denisova', 'vindija' ],
         [ 'altai', 'chag' ],
         [ 'altai', 'vindija' ],
         [ 'chag', 'vindija' ],
         [ 'denisova' ],
         [ 'altai' ],
         [ 'chag' ],
         [ 'vindija' ]]
    
    in_common_above0 = [0] * len(combinations)
    in_common_equal = [0] * len(combinations)
    in_common_arc = [0] * len(combinations)
    
    for combo in combinations:
        arc_name = combo[0]
        common_above0 = 0
        common_equal = 0
        common_arc = 0
        for window in range(len(archaics[arc_name])):
            val = archaics[arc_name][window]
            if all(archaics[arc][window] > 0 for arc in combo):
                common_above0 += 1
                if all(archaics[arc][window] == val for arc in combo):
                    common_equal += 1
            if all(archaics[arc][window] == 'Archaic' for arc in combo):
                common_arc += 1
        idx = combinations.index(combo)
        in_common_above0[idx] = common_above0
        in_common_equal[idx] = common_equal
        in_common_arc[idx] = common_arc
        
        # print(f'combo: {combo}, fraction of windows above 0: {common_above0 / len(kmer_counts)}')
        # print(f'combo: {combo}, fraction of windows equal: {common_equal / len(kmer_counts)}')
        # print(f'combo: {combo}, fraction of windows equal: {common_arc / len(kmer_counts)}')
            
    plot_data = from_memberships(
        [[ 'denisova', 'altai', 'chag', 'vindija' ],
         [ 'denisova', 'altai', 'chag' ],
         [ 'denisova', 'altai', 'vindija' ],
         [ 'denisova', 'chag', 'vindija' ],
         [ 'altai', 'chag', 'vindija' ],
         [ 'denisova', 'altai' ],
         [ 'denisova', 'chag' ],
         [ 'denisova', 'vindija' ],
         [ 'altai', 'chag' ],
         [ 'altai', 'vindija' ],
         [ 'chag', 'vindija' ],
         [ 'denisova' ],
         [ 'altai' ],
         [ 'chag' ],
         [ 'vindija' ]],
        data=in_common_equal
        )
    plot(plot_data, show_counts='%d')
    plt.title('Kmer counts equal')
    output = filename.replace('probs_and_path.tsv', 'kmers_equal.png')
    plt.savefig(output, dpi=300)


if __name__ == "__main__":
    main()
