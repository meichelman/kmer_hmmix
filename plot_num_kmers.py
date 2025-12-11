#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


def get_counts(file):
    kmer_counts = []
    
    with open(file) as f:
        
        for line in f:
            
            contig, start, end, count, introgressed = line.strip().split('\t')
            kmer_counts.append(int(count))

    return np.array(kmer_counts)


def plot(data, path):
    assemblies = list(data.keys())
    archaics = list(next(iter(data.values())).keys())
    num_assemblies = len(assemblies)
    num_archaics = len(archaics)

    # # Compute total k-mer counts
    # kmer_counts = np.zeros((num_assemblies, num_archaics), dtype=int)
    # for i, asm in enumerate(assemblies):
    #     for j, arc in enumerate(archaics):
    #         kmer_counts[i, j] = np.sum(data[asm][arc])  # sum of values in list

    # # Plot
    # fig, ax = plt.subplots(figsize=(10, 6))

    # # X positions for assemblies
    # x = np.arange(num_assemblies)

    # # Offset for archaics within each assembly column
    # width = 0.15
    # for j, arc in enumerate(archaics):
    #     ax.scatter(x + (j - num_archaics/2)*width, kmer_counts[:, j],
    #                color=plt.cm.tab10(j), label=arc, s=80, marker='o')

    # # Formatting
    # ax.set_xticks(x)
    # ax.set_xticklabels(assemblies, rotation=45, ha="right")
    # ax.set_ylabel("Total k-mer counts")
    # ax.set_title("K-mer counts per Assembly and Archaic")
    # ax.legend(title="Archaics")

    # plt.tight_layout()
    # plt.savefig(path + 'kmer_counts_scatter_plot.png', dpi=300)
    


    thresholds = [0, 10, 20, 30, 40]
    fractions = np.zeros((num_assemblies, num_archaics, len(thresholds)))

    for i, asm in enumerate(assemblies):
        for j, arc in enumerate(archaics):
            counts = data[asm][arc]
            total_windows = len(counts)
            for k, t in enumerate(thresholds):
                fractions[i, j, k] = np.sum(counts > t) / total_windows

    # --- NEW: Plot fractions ---
    fig, ax = plt.subplots(figsize=(10, 6))

    # Define AFR vs non-AFR mapping
    afr_assemblies = {"HG01891_h1", "HG02055_h1"}  # mark AFR explicitly
    color_map = {True: "red", False: "blue"}       # AFR = red, non-AFR = blue

    # Define line styles for individuals
    line_styles = ["-", "--", "-.", ":"]
    style_map = {asm: line_styles[i % len(line_styles)] for i, asm in enumerate(assemblies)}

    # Plot
    for j, arc in enumerate(archaics):
        for i, asm in enumerate(assemblies):
            is_afr = asm in afr_assemblies
            ax.plot(thresholds, fractions[i, j, :],
                    marker='o',
                    color=color_map[is_afr],
                    linestyle=style_map[asm],
                    label=asm if j == 0 else "")

    ax.set_xticks(thresholds)
    ax.set_xlabel("Threshold (k-mer count > t)")
    ax.set_ylabel("Fraction of windows")
    ax.set_title("Fraction of windows above thresholds per Assembly and Archaic")

    # Legend will show each assembly once
    ax.legend(title="Assembly (AFR=red, non-AFR=blue)", bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(path + 'fraction_windows_thresholds.png', dpi=300)


def main():
    
    root = '9sub_runs/'
    assemblies = [ 'HG01891_h1', 'HG02055_h1', 'HG002_h1', 'HG00438_h1', 'HG005_h1' ]
    # AFR, AFR, EUR, EAS, EAS
    archaics = [ 'denisova', 'altai', 'chag', 'vindija' ]
    data = {asm: {arc: [] for arc in archaics} for asm in assemblies}
    
    for asm in assemblies:    
        for arc in data[asm].keys():
            filename = 'windows_across_genome_with_zero_and_nonzero_matching_kmers_and_including_introgressed_and_no_introgressed_regions.bed'
            filepath = root + asm + '/' + arc + '/' + filename
            kmer_counts = get_counts(filepath)
            print(f'Assembly: {asm}, Archaic: {arc}, max count in a window: {np.max(kmer_counts)}')
            # data[asm][arc] = kmer_counts
    
    # plot(data, root)
    

if __name__ == "__main__":
    main()
