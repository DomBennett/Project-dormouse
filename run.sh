#!/bin/bash
# Dom Bennett
# Take BLAST results, output tree
# Requires: mafft-qinsi, raxml, python2.7 and biopython
# run by typing: `sh run.sh >& log.txt &`

echo Started: $(date)

# EXECUTABLE PATHS
MAFFT=/home/djb208/bin/mafft-qinsi
RAXML=$(which raxml)
PYTHON=$(which python)

# FILES AND FOLDERS
WD=$(pwd)
ALIGN=${WD}/2_alignments/supermatrix.phy
TREEOUT=${WD}/3_trees
PART=${WD}/2_alignments/partitions.txt

# SORT SEQUENCES
echo '\n1. Sorting sequences ....'
$PYTHON python_scripts/sort_sequences.py
echo 'Complete\n'

# RUN ALIGNMENTS
echo '2. Running MAFFT for ....'
if [ ! -d 2_alignments ]; then
  mkdir 2_alignments
fi
for SF in `ls 1_sequences | grep \.fasta`
do
  echo ".... $SF"
  $MAFFT 1_sequences/$SF > 2_alignments/${SF%.fasta}_alignment.fasta 2> 2_alignments/${SF%.fasta}_log.txt
done
echo 'Complete\n'

# COMBINE ALIGNMENTS
echo '3. Combining alignments ....'
$PYTHON python_scripts/combine.py
echo 'Complete\n'

# RUN RAXML
if [ ! -d 3_trees ]; then
  mkdir 3_trees
fi
echo '4. Running RAXML with bootstrap'
$RAXML -n S1 -m GTRCAT -p $RANDOM -# 100 -s $ALIGN -T 2 -q $PART  -o outgroup -w $TREEOUT
$RAXML -n S2 -m GTRCAT -p $RANDOM -b $RANDOM -# 100 -s $ALIGN -T 2 -q $PART -o outgroup -w $TREEOUT
$RAXML -n S3 -m GTRCAT -p $RANDOM -f b -t ${TREEOUT}/RAxML_bestTree.S1 -z ${TREEOUT}/RAxML_bootstrap.S2 -w $TREEOUT
echo 'Complete\n'

echo Finished: $(date)
