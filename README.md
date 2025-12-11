## Matt Eichelman 2nd Rotation Materials
This page contains the materials for using the HMM implemented during the rotation.

In order to use the HMM for training and/or decoding the following is required:
- Output from ARCkmerFinder
- Directory named after the assembly run through ARCkmerFinder (i.e. HG002_h1)
    - Subdirectories named 'altai', 'denisova', 'chag', 'vindija'
- Conda environment with the packages listed in env.yml

Previous runs can be found in the 9sub_runs directory

Walkthrough of the pipeline:
1. Run ARCkmerFinder with the assembly and save the file named 'windows_across_genome_with_zero_and_nonzero_matching_kmers_and_including_introgressed_and_no_introgressed_regions.bed'
2. Create a directory named after the assembly of interest (i.e. HG002_h1)
3. Create subdirectories altai, denisova, chag, vindija
4. Create and activate the conda environment
5. Run the train_decode_analyze_pipeline.py script
    - './train_decode_analyze_pipeline.py [assembly name]'
    - Example: './train_decode_analyze_pipeline.py HG002_h1'

This script creates the mutation rate file, trained parameters file, and decoding file along with a couple plots I used for some testing.

The output files can be found in the subdirectory mutrate_eq1MbBin.

For running individual parts of the pipeline like generating the mutation rates, training, and decoding there is an example in main.py (commented out in the script, and shown below). Additionally, the code in train_decode_analyze_pipeline.py makes individual calls so the lines can be modified to use on the command line.

For making individual calls (must be done in this order):
1. Generate mutation rates
    - python3 main.py mutrate HG01891_h1/vidija/HG01891_h1_counts.bed HG01891_h1/vidija/mutrate_eq1MbBin/HG01891_h1_mutrates.txt
2. Training
    - python3 main.py train none HG01891_h1/vidija/HG01891_h1_counts.bed HG01891_h1/vidija/mutrate_eq1MbBin/HG01891_h1_mutrates.txt HG01891_h1/vidija/mutrate_eq1MbBin/trained.json
    - The 'none' argument was for testing and MUST be included
3. Decoding
    - python3 main.py decode HG01891_h1/vidija/mutrate_eq1MbBin/trained.json HG01891_h1/vidija/HG01891_h1_counts.bed HG01891_h1/vidija/mutrate_eq1MbBin/HG01891_h1_mutrates.txt HG01891_h1/vidija/mutrate_eq1MbBin/probs_and_path.tsv