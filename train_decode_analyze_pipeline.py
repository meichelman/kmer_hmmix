#!/usr/bin/env python

import subprocess
import sys
import os

"""
Steps to do before running this script:
1. Activate conda env: "conda activate hmmix"
2. Ensure all directories and files exist (assembly, archaics, counts.bed)
"""


def main():
    args = sys.argv
    if len(args) != 2:
        sys.exit('\n\nWrong input\n\n')
        
    assembly_path = args[1]              # full path, e.g. "9sub_runs/HG002_h1"
    assembly_name = os.path.basename(assembly_path)  # just "HG002_h1"
    
    archaics = ['altai', 'denisova', 'chag', 'vindija']
    for arc_dir_name in os.listdir(assembly_path):
        if arc_dir_name not in archaics:
            continue
        arc_dir = os.path.join(assembly_path, arc_dir_name)
        if os.path.isdir(arc_dir):
            
            counts_bed = None
            bed_filename = 'windows_across_genome_with_zero_and_nonzero_matching_kmers_and_including_introgressed_and_no_introgressed_regions.bed'
            for file in os.listdir(arc_dir):
                if file == bed_filename:
                    counts_bed = os.path.join(arc_dir, file)
            if counts_bed is None:
                sys.exit(f'\n\nERROR! counts.bed file not found in {arc_dir}\n\n')

            # Mutation rate step
            mut_out = os.path.join(arc_dir, 'mutrate_eq1MbBin', f'{assembly_name}_mutrates.bed')
            mutrate_command = [
                "python",
                "main.py",
                "mutrate",
                counts_bed,
                mut_out
            ]
            # print(f'About to execute: {mutrate_command}')
            # subprocess.run(mutrate_command, check=True)

            # Training step
            train_out = os.path.join(arc_dir, 'mutrate_eq1MbBin', 'trained.json')
            train_command = [
                "python",
                "main.py",
                "train",
                "none",
                counts_bed,
                mut_out,
                train_out
            ]
            # print(f'About to execute: {train_command}')
            # subprocess.run(train_command, check=True)

            # Decoding step
            decode_out = os.path.join(arc_dir, 'mutrate_eq1MbBin', 'probs_and_path.tsv')
            decode_command = [
                "python",
                "main.py",
                "decode",
                train_out,
                counts_bed,
                mut_out,
                decode_out
            ]
            # print(f'About to execute: {decode_command}')
            # subprocess.run(decode_command, check=True)
            
            # Analysis step
            plot_kmers_vs_probs_command = [
                './plot_kmers_vs_probs.py',
                decode_out
            ]
            # print(f'About to execute: {plot_kmers_vs_probs_command}')
            # subprocess.run(plot_kmers_vs_probs_command, check=True)
            
            plot_consecutive_arc_states_command = [
                './plot_consecutive_arc_states.py',
                decode_out
            ]
            print(f'About to execute: {plot_consecutive_arc_states_command}')
            subprocess.run(plot_consecutive_arc_states_command, check=True)
            
    
    
if __name__ == "__main__":
    main()
