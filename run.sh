#!/bin/bash
# Dom Bennett
# Take BLAST results, output tree
# Requires: mafft-qinsi, raxml, python2.7 and biopython
# run by typing: `sh run.sh >& log.txt &`

# EXECUTABLE PATHS
MAFFT=/home/djb208/bin/mafft-qinsi
RAXML=/usr/biosoft/bin/raxml
PYTHON=/usr/local/python/2.7.5/bin/python

# FILES AND FOLDERS
WD=$pwd
ALIGN=$WD/2_alignments/supermatrix.phy
TREEOUT=$WD/3_trees

# SORT SEQUENCES
$PYTHON python_scripts/sort_sequences.py

# RUN ALIGNMENTS
if [ ! -d 2_alignments ]; then
  mkdir 2_alignments
fi
echo 'Running MAFFT'
echo '.... p1'
$MAFFT 1_sequences/p1.fasta > 2_alignments/p1_alignment.fasta
echo '.... p2'
$MAFFT 1_sequences/p2.fasta > 2_alignments/p2_alignment.fasta
echo '.... p3'
$MAFFT 1_sequences/p3.fasta > 2_alignments/p3_alignment.fasta
echo 'Complete'

# COMBINE ALIGNMENTS
$PYTHON python_scripts/combine.py

# RUN RAXML
if [ ! -d 3_trees ]; then
  mkdir 3_trees
fi
echo 'Running RAXML'
$RAXML -s $ALIGN -n tree -m GTRCAT -o outgroup -w $TREEOUT -T 2 -p $RANDOM
echo 'Complete'
