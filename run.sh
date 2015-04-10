#!/bin/bash
# Dom Bennett
# Take BLAST results, output tree
# Requires: mafft-qinsi, python2.7 and biopython

# EXECUTABLE PATHS
$MAFFT='~/bin/mafft-qinsi'
$RAXML='raxml'
$PYTHON='/bin/usr/python2.7'

# FILES AND FOLDERS
$ALIGN = pwd/2_alignments/supermatrix.phy
$TREEOUT = pwd/3_trees

# SORT SEQUENCES
$PYTHON python_scripts/sort_sequences.py

# RUN ALIGNMENTS
if [ ! -d 2_alignments ]; then
  mkdir 2_alignments
fi
echo 'Running MAFFT'
echo '.... p1'
$MAFFT p1.fasta > 2_alignments/p1_alignment.fasta
echo '.... p2'
$MAFFT p2.fasta > 2_alignments/p2_alignment.fasta
echo '.... p3'
$MAFFT p3.fasta > 2_alignments/p3_alignment.fasta
echo 'Complete'

# COMBINE ALIGNMENTS
$PYTHON python_scripts/combine.py

# RUN RAXML
if [ ! -d 3_trees ]; then
  mkdir 3_trees
fi
echo 'Running RAXML'
$RAXML -s $ALIGN -n tree -m GTRCAT -o outgroup -w $TREEOUT
echo 'Complete'
