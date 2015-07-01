#!/bin/bash
# Dom Bennett
# Take BLAST results, output tree
# Requires: mafft-qinsi, raxml, python2.7 and biopython
# run by typing: `sh run_raxml.sh >& log_raxml.txt &`

echo Started: $(date)

# EXECUTABLE PATHS
RAXML=$(which raxml)

# FILES AND FOLDERS
WD=$(pwd)
ALIGN=${WD}/2_alignments/supermatrix.phy
TREEOUT=${WD}/3_trees/
PART=${WD}/2_alignments/partitions.txt

# add -o outgroup and -q $PART if need be
echo '4. Running RAXML with bootstrap'
$RAXML -n S1 -m GTRCAT -p $RANDOM -# 100 -s $ALIGN -w $TREEOUT
$RAXML -n S2 -m GTRCAT -p $RANDOM -b $RANDOM -# 100 -s $ALIGN -w $TREEOUT
$RAXML -n S3 -m GTRCAT -p $RANDOM -f b -t ${TREEOUT}RAxML_bestTree.S1 -z ${TREEOUT}RAxML_bootstrap.S2 -w $TREEOUT
echo 'Complete\n'

echo Finished: $(date)
