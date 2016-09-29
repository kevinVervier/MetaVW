#Copyright (c) by respective owners including bioMerieux, and
#individual contributors Kevin Vervier and Pierre Mahe. All #rights reserved.  Released under a BSD (revised)
#license as described in the file LICENSE.

# specify path to fasta2vw
fasta2vw=../../../tools/fasta2vw


# specify reference database
#--------------------------#
DB=small
#DB=large


# specify parameters #
#--------------------#
NBATCHES=10
K=12
modelDir=../../2-build-models/output/train_$DB-db
model=$modelDir/vw-modelvw-model_batch-$NBATCHES.model
dico=$modelDir/vw-dico.txt


# make predictions fragments #
#----------------------------#
fasta=../../1-generate-test-datasets/output/test.fragments.fasta
prefix=../output/test.fragments.$DB-db
# get vw predictions
$fasta2vw -i $fasta -k $K | vw -t -i $model -p $prefix.preds.vw
# convert vw class to taxid
Rscript vw-class-to-taxid.R $prefix.preds.vw $dico $prefix.preds.taxid


# make predictions - reads mutation #
#------------------------------------#
for freq in 1 2 3 4 5 6 7 8 9 10 10.7
do
	fasta=../../1-generate-test-datasets/output/test.mutation-$freq-reads.fastq
	prefix=../output/test.mutation-$freq-reads.$DB-db
	# get vw predictions
	$fasta2vw -i $fasta -k $K | vw -t -i $model -p $prefix.preds.vw
	# convert vw class to taxid
	Rscript vw-class-to-taxid.R $prefix.preds.vw $dico $prefix.preds.taxid
done


# make prediction - reads homo #
#------------------------------#
for homo in Balzer Richter Margulies
do
	fasta=../../1-generate-test-datasets/output/test.homo-$homo-reads.fastq
	prefix=../output/test.homo-$homo-reads.$DB-db
	# get vw predictions
	$fasta2vw -i $fasta -k $K | vw -t -i $model -p $prefix.preds.vw
	# convert vw class to taxid
	Rscript vw-class-to-taxid.R $prefix.preds.vw $dico $prefix.preds.taxid
done


