#Copyright (c) by respective owners including bioMerieux, and
#individual contributors Kevin Vervier and Pierre Mahe. All #rights reserved.  Released under a BSD (revised)
#license as described in the file LICENSE.


# specify reference database
#--------------------------#
DB=small

# specify parameters #
#--------------------#
NBATCHES=10
K=10
modelDir=../../2-build-models/output/train_$DB-db
model=$modelDir/vw-model_batch-$NBATCHES.model
dico=$modelDir/vw-dico.txt

# make predictions fragments #
#----------------------------#
L=200
fasta=../../../data/small/test-dataset/fragments_L-$L.fasta
prefix=../output/preds.fragments_L-$L.from-vw
# get vw predictions
fasta2vw -i $fasta -k $K | vw -t -i $model -p $prefix.preds.vw -r $prefix.scores

# convert vw class to taxid
Rscript vw-class-to-taxid.R $prefix.preds.vw $dico $prefix.preds.taxid


