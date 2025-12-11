import numpy as np
from collections import defaultdict

from helper_functions import sortby, Make_folder_if_not_exists

def make_mutation_rate(countfile, outfile, window_size):

    kmers_counts_window = defaultdict(lambda: defaultdict(int))
    contig_lengths = defaultdict(int)
    with open(countfile) as data:
        assembly = data.readline().split('#')[0]
        data.seek(0)
        for line in data:
            if line.startswith(assembly):
                # testing with kmer_hmmix/HG01891_h1/vidija/HG01891_h1_counts.bed
                if len(line.split('\t')) == 5:
                    contig, window_start, window_end, count, _ = line.split('\t')
                else:
                    contig, window_start, window_end, count = line.split('\t')
                window = int(window_start) - int(window_start) % window_size
                kmers_counts_window[contig][window] += int(count)
                contig_lengths[contig] = int(window_end)
    
    kmers = []
    genome_positions = []
    for contig in sorted(kmers_counts_window, key=sortby):
        lastwindow = max(kmers_counts_window[contig]) + window_size
        lastwindow = int(lastwindow)
        # print(f'Contig: {contig}, end: {contig_lengths[contig]}')

        for window in range(0, lastwindow, window_size):
            kmers.append(kmers_counts_window[contig][window])
            genome_positions.append([contig, window, window + window_size])

    kmers = np.array(kmers)
    assembly_length = sum(contig_lengths.values())
    # print(f'Assembly length: {assembly_length}')
    # Assembly length: 3043232268
    genome_mean = np.sum(kmers) / np.sum(assembly_length)
    # print(f'Genome mean: {genome_mean}')
    # Genome mean: 0.0007280896773141077

    Make_folder_if_not_exists(outfile)
    with open(outfile, 'w') as out:
        print('contig', 'start', 'end', 'mutationrate', sep = '\t', file = out)
        for genome_pos, mut in zip(genome_positions, kmers):
            contig, start, end = genome_pos
            if mut == 0:
                ratio = 0
            else:
                ratio = round(mut / window_size / genome_mean, 5)
                # 10 / 2000 / 0.0005 

            
            print(contig, start, end, ratio, sep = '\t', file = out)
            # print(contig, start, end, np.random.choice(range(1, 41)) , sep = '\t', file = out)
            # print(contig, start, end, 1, sep = '\t', file = out)
            # print(contig, start, end, genome_mean * 2000, sep = '\t', file = out)
