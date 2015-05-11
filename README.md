## Phylogenetic analysis of dormouse tapeworm samples

Pipeline for determining species of tapeworm discovered in Hazel dormouse (*Muscardinus avellanarius*) using 18S and 28S PCR amplified sequences.

### Stages

1. NCBI BLAST amplified sequences against Cyclophellyidea species (txid6201) [MANUALLY]
2. NCBI BLAST amplified sequences against Caryophyllidea species (txid108240) [MANUALLY]
3. Align with MAFFT-QINISI each amplified sequence against BLAST results
4. Combine into single supermatrix
5. Use RAxML to create tree with Caryophyllidea outgroup and partition by gene

The scripts provided here perform stages 3-4. Original data not provided.

### Requirements

* Python
* Biopython
* MAFFT
* RAxML
* UNIX

### Running

```{bash}
run.sh >& log.txt
```

### Dir structure
```
-- run.sh
-- 1_sequences/
-- 2_alignments/
-- 3_trees/
-- python_scripts/
  -- sort_names.py
    -- [converts the names of each in the blast results to species names for tree]
  -- combine.py
    -- [combines the results of each alignment BLAST results intoa single supermatrix]
-- data/
  -- samples/
    -- [PCR amplified sequences from dormouse sample]
  -- BLAST_results/
    -- [Top 100 BLAST results for each amplified sequence constrained to the Cyclophyllidea, and then Caryophyllidea (outgroup). We used megablast with default parameters.]
-- misc/
    -- all_cyclophyllidea_18s_sequences.fasta
      -- [Downloaded all 18S sequences for Cyclophyllidea (tapeworms) NCBI nucleotide search term: "txid6201[Organism:exp] AND 18S NOT predicted[TI] NOT shotgun[TI] NOT scaffold[TI] NOT assembly[TI] NOT unverified[TI]"]
    -- primer_pairs.fasta
      -- [the primer pairs for 18S that Gaby sequenced from the Dormouse sample]

```

### Run history
BLAST results were downloaded from NCBI manually from 04-05/2015

run.sh was run on Imperial College London's bioinformatics server: codon.bioinformatics.ic.ac.uk (05/2015)

### Author
Dom Bennett
