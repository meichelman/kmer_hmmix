import numpy as np
import matplotlib.pyplot as plt

# --- Load counts from BED file ---
bed = 'HG01891_h1/denisova/HG01891_h1_counts.bed'
counts = []
with open(bed) as f:
    for line in f:
        id, start, end, count, _ = line.strip().split('\t')
        # if count == '1173':
        #     print(line)
        counts.append(int(count))  # convert to integer

counts = np.array(counts)
# print(f'Counts array length: {len(counts)}')
# filtered_counts = counts[counts > 50]
# print(len(filtered_counts) / len(counts) * 100)

# --- Remove zeros ---
counts_no_zeros = counts[counts != 0]

# --- Exclude outliers using IQR ---
Q1 = np.percentile(counts_no_zeros, 25)
Q3 = np.percentile(counts_no_zeros, 75)
IQR = Q3 - Q1
lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR

# filtered_counts = counts_no_zeros[(counts_no_zeros >= lower_bound) & (counts_no_zeros <= upper_bound)]
filtered_counts = counts_no_zeros

# --- Compute summary statistics (zeros + outliers removed) ---
summary = {
    "count": len(filtered_counts),
    "mean": np.mean(filtered_counts),
    "median": np.median(filtered_counts),
    "std_dev": np.std(filtered_counts, ddof=1),
    "min": np.min(filtered_counts),
    "max": np.max(filtered_counts),
    "Q1 (25th percentile)": np.percentile(filtered_counts, 25),
    "Q3 (75th percentile)": np.percentile(filtered_counts, 75),
    "IQR": np.percentile(filtered_counts, 75) - np.percentile(filtered_counts, 25)
}

print("Summary statistics (zeros removed, outliers removed):")
for k, v in summary.items():
    print(f"{k}: {v}")

# --- Violin plot (zeros + outliers removed) ---
# plt.figure(figsize=(8,6))
# plt.violinplot(filtered_counts, showmeans=True, showmedians=True)

# plt.ylabel("Counts", fontsize=14)
# plt.title("Distribution of Counts (Zeros & Outliers Removed)", fontsize=16, weight='bold')
# plt.xticks([1], ["Counts"])
# plt.grid(axis='y', linestyle='--', alpha=0.7)
# plt.tight_layout()
# plt.show()
