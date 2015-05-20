#!/bin/bash
# Dom Bennett
# Take BLAST results, output tree
# Requires: mafft-qinsi, raxml, python2.7 and biopython
# run by typing: `sh run.sh >& log.txt &`

# EXECUTABLE PATHS
MAFFT=/home/djb208/bin/mafft-qinsi
RAXML=$(which raxml)
PYTHON=$(which python)

# FILES AND FOLDERS
WD=$(pwd)
ALIGN=${WD}/2_alignments/supermatrix.phy
TREEOUT=${WD}/3_trees
PART=${WD}/2_alignments/partitions.txt

echo $WD
echo $PART
echo $ALIGN
echo $TREEOUT

# SORT SEQUENCES
$PYTHON python_scripts/sort_sequences.py

# RUN ALIGNMENTS
if [ ! -d 2_alignments ]; then
  mkdir 2_alignments
fi
for seq_file in `ls 1_sequences | grep \.fasta`
do
  echo "Running MAFFT for $seq_file ...."
  $MAFFT 1_sequences/$seq_file > 2_alignments/${seq_file}_alignment.fasta 2> 2_alignments/${seq_file}_alignment_log.txt
done
echo 'Complete'

# # COMBINE ALIGNMENTS
$PYTHON python_scripts/combine.py

# RUN RAXML
if [ ! -d 3_trees ]; then
  mkdir 3_trees
fi
echo 'Running RAXML with bootstrap'
$RAXML -n S1 -m GTRCAT -p $RANDOM -# 100 -s $ALIGN -T 2 -q $PART  -o outgroup -w $TREEOUT
$RAXML -n S2 -m GTRCAT -p $RANDOM -b $RANDOM -# 100 -s $ALIGN -T 2 -q $PART -o outgroup -w $TREEOUT
$RAXML -n S3 -m GTRCAT -p $RANDOM -f b -t ${TREEOUT}/RAxML_bestTree.S1 -z ${TREEOUT}/RAxML_bootstrap.S2 -w $TREEOUT
echo 'Complete'
