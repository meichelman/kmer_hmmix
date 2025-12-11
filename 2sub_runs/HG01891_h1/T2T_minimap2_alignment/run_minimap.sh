#!/bin/bash -l


#SBATCH --time=24:00:00
#SBATCH --mem=120g
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eiche118@umn.edu

set -v

module load samtools/1.14

szMinimapExe=/projects/standard/hsiehph/shared/software/packages/minimap2/minimap2
szReference=/projects/standard/hsiehph/shared/hsiehph_shared/assemblies/CHM13/T2T/v2.0/ucsc_fixed_chr_order/hs1.fa
szQuery=/projects/standard/hsiehph/gordo893/Matt/T2T_minimap2_alignment/HG01891#1#JAGYVO010000020.1:12285000-12385000.fa
szOutput=alignments.paf
# --eqx and -c are necessary to use rustybam liftover
$szMinimapExe -t 10 --eqx -c --cs -x asm5 --secondary=no -s 25000 -K 8G $szReference $szQuery >$szOutput
szOutput2=alignments.bed
cat $szOutput | awk '{print $6"\t"$8"\t"$9"\t"$1"\t"$3"\t"$4"\t"$5 }' >$szOutput2
