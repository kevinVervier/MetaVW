
#Copyright (c) by respective owners including bioMerieux, and
#individual contributors Kevin Vervier and Pierre Mahe. All #rights reserved.  Released under a BSD (revised)
#license as described in the file LICENSE.

# specify path to drawfrag
drawfrag=../../../tools/drawfrag

# specify input  data
fasta=../../../data/test-dataset/test_db.fasta
taxids=../../../data/test-dataset/test_db.species-level.taxid

# specify fragments parameters (length and cover ~ number)
L=200
COVERAGE=1

#############################
###### SMALL SCALE TEST #####
#COVERAGE=0.05
#############################

# set seed (for reproducibility)
SEED=42

# draw fragments
$drawfrag -i $fasta -t $taxids -l $L -c $COVERAGE -o ../output/test.fragments.fasta -g ../output/test.fragments.gi2taxid -s $SEED

# extract taxids
cut -f 2 ../output/test.fragments.gi2taxid > ../output/test.fragments.taxid


