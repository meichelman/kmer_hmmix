#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
# import sys


# kmer_file, GC_file = sys.argv[1], sys.argv[2]
# kmer_counts = []
# GC_content = []
# starts = []
# ends = []

# Parse the file

# with open(kmer_file) as infile1, open(GC_file) as infile2:
#     for line in infile1:
#         if line.startswith('start'):
#             continue
#         # start, end, count, mutrate, human_p, arc_p, post_state, viterbi_state = line.strip().split('\t')
#         contig, start, end, count = line.strip().split('\t')
#         kmer_counts.append(int(count))
#         starts.append(int(start))
#         ends.append(int(end))
        
#     for line in infile2:
#         contig, start, end, gc = line.strip().split('\t')
#         GC_content.append(float(gc))

# kmer_counts = np.array(kmer_counts)
# # print(len(kmer_counts))
# # print(len(GC_content))
# # GC_content_contig1 = np.array(GC_content[:55638])  # Ensure same length


# def view_and_zoom():
#     starts = np.array(starts)
#     ends = np.array(ends)

#     # Convert bp to Mb for x-axis
#     start_mb_array = starts / 1e6
#     end_mb_array = ends / 1e6

#     # ---- Zoom window definition ----
#     start_mb, end_mb = 10, 15  # Mb

#     # Interval overlap mask: window overlaps region if:
#     #   interval_end ≥ zoom_start  AND  interval_start ≤ zoom_end
#     mask = (end_mb_array >= start_mb) & (start_mb_array <= end_mb)

#     # Apply mask
#     x_mb_zoom = start_mb_array[mask]
#     kmer_counts_zoom = kmer_counts[mask]
#     GC_content_zoom = GC_content[mask]

#     # ---- Plot ----
#     fig, ax1 = plt.subplots(figsize=(12, 5))

#     # Left y-axis: kmer count
#     ax1.set_xlabel("Genomic position (Mb)")
#     ax1.set_ylabel("kmer count", color='tab:blue')
#     ax1.plot(x_mb_zoom, kmer_counts_zoom, color='tab:blue', linewidth=1)
#     ax1.tick_params(axis='y', labelcolor='tab:blue')

#     # Right y-axis: archaic probability (faded)
#     ax2 = ax1.twinx()
#     ax2.set_ylabel("GC content", color='tab:red')
#     ax2.plot(x_mb_zoom, GC_content_zoom, color='tab:red', linewidth=1, alpha=0.4)
#     ax2.tick_params(axis='y', labelcolor='tab:red')

#     plt.tight_layout()
#     output = kmer_file.replace('probs_and_path.tsv', 'kmers_vs_GC.png')
#     plt.savefig(output, dpi=300)


def get_counts(counts, GC):
    kmer_counts = []
    GC_content = []
    
    with open(counts) as f1, open(GC) as f2:

        for line1, line2 in zip(f1, f2):
            contig1, start1, end1, count, introgressed = line1.strip().split('\t')
            contig2, start2, end2, gc = line2.strip().split('\t')

            count = int(count)
            gc = float(gc)

            if count > 0:
                kmer_counts.append(count)
                GC_content.append(gc)

    return np.array(kmer_counts), np.array(GC_content)


def plot_correl(data, path):
    assemblies = list(data.keys())
    archaics = list(next(iter(data.values())).keys())
    num_assemblies = len(assemblies)
    num_archaics = len(archaics)

    # Build correlation matrix from data dict
    correl_matrix = np.zeros((num_assemblies, num_archaics))
    for i, asm in enumerate(assemblies):
        for j, arc in enumerate(archaics):
            correl_matrix[i, j] = data[asm][arc]

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(num_assemblies)
    width = 0.15

    for j, arc in enumerate(archaics):
        ax.scatter(x + (j - num_archaics/2)*width, correl_matrix[:, j],
                   color=plt.cm.tab10(j), label=arc, s=80, marker='o')

    ax.set_xticks(x)
    ax.set_xticklabels(assemblies, rotation=45, ha="right")
    ax.set_ylabel("Correlation coefficient")
    ax.set_title("Correlation between k-mer counts and GC content")
    ax.legend(title="Archaics")
    
    plt.tight_layout()
    plt.savefig(path + 'kmers_GC_correl.png', dpi=300)


def plot_points(GC_content, kmer_counts, path):
    plt.plot(kmer_counts, GC_content, 'o', markersize=1, alpha=0.5)
    plt.xlabel("GC content")
    plt.ylabel("kmer counts")
    plt.title(f"HG005_h1, vindija database")
    plt.tight_layout()
    plt.savefig(path + '/kmer_GC_vs_counts_HG005_h1_vindija.png', dpi=300)

    
def main():
    # view_and_zoom()
    # correlation_test()
    
    root = '9sub_runs/'
    assemblies = [ 'HG01891_h1', 'HG02055_h1', 'HG002_h1', 'HG00438_h1', 'HG005_h1' ]
    # AFR, AFR, EUR, EAS, EAS
    archaics = [ 'denisova', 'altai', 'chag', 'vindija' ]
    data = {asm: {arc: [] for arc in archaics} for asm in assemblies}
    
    for asm in assemblies:
        GC_path = root + asm + f'/{asm}_GC_content.tsv'
        for arc in data[asm].keys():
            filename = 'windows_across_genome_with_zero_and_nonzero_matching_kmers_and_including_introgressed_and_no_introgressed_regions.bed'
            filepath = root + asm + '/' + arc + '/' + filename
            kmer_counts, GC_content = get_counts(filepath, GC_path)
            correl = np.corrcoef(GC_content, kmer_counts)[0, 1]
            data[asm][arc] = round(correl, 4)
            
    path = root + asm + '/' + arc
    plot_points(kmer_counts, GC_content, path)
    plot_correl(data, root)
    

if __name__ == "__main__":
    main()
