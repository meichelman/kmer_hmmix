
# obtain assemblies
...

# count kmers and save to a db
meryl count k=21 /projects/standard/hsiehph/shared/hsiehph_shared/assemblies/HPRCasm/HG01891_h1.fa.gz output HG01891_h1.meryl
meryl -t 1 -Q print arc_plus_afr.meryl | head

# count kmers in other assemblies abd combine all into one db
meryl \ 
    union-sum output=union.meryl \
    [count /path/to/fna.gz output=assem2.meryl] \
    [count /path/to/fna.gz output=assem3.meryl] \
    assem1.meryl

# histogram of kmers
meryl -Q histogram assem1.meryl

