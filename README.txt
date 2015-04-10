# Phylogenetic analysis of dormouse tapeworm samples
# Dom Bennett
# 27/03/2015

# RUN HISTORY
BLAST results were downloaded from NCBI manually.
Sequences in sequences/ were compiled manually.
run.sh was run on codon.bioinformatics.ic.ac.uk ##/04/2015

# DIR STRUCTURE
-run.sh
    *run the analysis from alignment of sequences to trees
—sequences/
    *Three sequence fasta files for each pair containing BLAST results and
     primer pair sequences
—results/
    *contains all the produced trees
-python_scripts/
    —sort_names.py
        *converts the names of each in the blast results to species names for tree
    -combine.py
        *combines the results of each primer pair's alignment BLAST results into
         a single supermatrix
-data/
    - p1.fasta
        * dormouse sample primer pair 1
    - p2.fasta
        * dormouse sample primer pair 2
    - p3.fasta
        * dormouse sample primer pair 3
    — BLAST_results/
        *Top 100 BLAST results for each pair constrained to the Cyclophyllidea
        *Top BLAST results for each pair constrained to the Caryophyllidea
         (outgroup)
        *We used megablast with default parameters.

    -misc/
        — all_cyclophyllidea_18s_sequences.fasta
            *Downloaded all 18S sequences for Cyclophyllidea (tapeworms)
             NCBI nucleotide search term: "txid6201[Organism:exp] AND 18S NOT
             predicted[TI] NOT shotgun[TI] NOT scaffold[TI] NOT assembly[TI] NOT
             unverified[TI]"
        — primer_pairs.fasta
            *the primer pairs that Gaby sequenced from Dormouse sample
