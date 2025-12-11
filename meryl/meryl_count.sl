#!/bin/bash -l
#SBATCH --job-name=meryl_count
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00 # hrs:min:sec
#SBATCH --partition=msibigmem,msilarge,msismall
#SBATCH --mail-type=END,FAIL # (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=eiche118@umn.edu
#SBATCH --output=/users/0/eiche118/out/meryl_count.out
#SBATCH --error=/users/0/eiche118/err/meryl_count.err


# Define global vars
READS_DIR="/projects/standard/hsiehph/shared/hsiehph_shared/assemblies/HPRCasm"
DB_DIR="/users/0/eiche118/kmer_hmm/afr_run_no_subtract/meryl_dbs"
IDS_FILE="${READS_DIR}/IDfiles/AFR.new.ids"
NTHREADS=4

config() {
    set -euo pipefail
}

main() {
    IDS=$(grep -oE 'HG[0-9]+|NA[0-9]+' "$IDS_FILE")
    for ID in $IDS; do
        for hap in h1 h2; do
            READS="${READS_DIR}/${ID}_${hap}.fa.gz"
            DB="${DB_DIR}/${ID}_${hap}.meryl"

            echo "Processing $READS -> $DB"
            meryl count k=21 memory=16 threads=$NTHREADS $READS output $DB
        done
    done
}

main
