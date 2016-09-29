#!/bin/bash

#Copyright (c) by respective owners including bioMerieux, and
#individual contributors Kevin Vervier and Pierre Mahe. All #rights reserved.  Released under a BSD (revised)
#license as described in the file LICENSE.

######################################################
# SPECIFY REFERENCE DATABASE TO USE (small or large) #
######################################################
DB=small
######################################################


# specify path to drawfrag and fasta2vw #
#---------------------------------------#
drawfrag=../../../tools/drawfrag
fasta2vw=../../../tools/fasta2vw


# specify input data #
#--------------------#
fasta=../../../data/train-dataset/train_$DB-db.fasta
taxids=../../../data/train-dataset/train_$DB-db.species-level.taxid
# extract number of labels
NLABELS=$(sort -u $taxids | wc -l)

# specify fragments parameters #
#------------------------------#
NBATCHES=10
COVERAGE=1
L=200


#############################
#### SMALL COVERAGE TEST ####
#############################
#COVERAGE=0.05
#NBATCHES=3
#############################


# specify kmer size #
#-------------------#
K=12

# specify generic VW options #
#----------------------------#
LAMBDA1=0
LAMBDA2=0
NPASSES=1
BITS=31


# specify output data #
#---------------------#
outputDir=../output/train_$DB-db
mkdir $outputDir
# define output "dictionary" : taxid <--> vw classes
dico=$outputDir/vw-dico.txt
# define model prefix
modelPrefix=$outputDir/vw-model


# specify temporary directory (to shuffle input files) #
#------------------------------------------------------#
TMPDIR=../output/TMP/
mkdir $TMPDIR

# loop on batches #
#-----------------#
SEED=42
for i in $(seq 1 ${NBATCHES})
do
	# modify random seed
        SEED=$(($SEED+$i))

	#  draw fragments
	fastaBatch=$outputDir/train.batch-$i.fasta
	gi2taxidBatch=$outputDir/train.batch-$i.gi2taxid
	taxidBatch=$outputDir/train.batch-$i.taxid
	$drawfrag -i $fasta -t $taxids -l $L -c $COVERAGE -o $fastaBatch -g $gi2taxidBatch -s $SEED
	# extract taxids
	cut -f 2 $gi2taxidBatch > $taxidBatch

	# learn model
	if [[ $i -eq 1 ]]; then
		# first iteration : no previous model to build upon 		
		$fasta2vw -i $fastaBatch -t $taxidBatch -k $K -d $dico | awk 'BEGIN{srand($SEED);} {printf "%06d %s\n", rand()*1000000, $0;}' | sort -n -T $TMPDIR | cut -c8- | vw --random_seed $SEED --passes $NPASSES -f ${modelPrefix}_batch-${i}.model.sr --oaa $NLABELS -b $BITS --l1 $LAMBDA1 --l2 $LAMBDA2 --save_resume --progress 100000
	else
		# futher iterations : save-resume mechanism
		$fasta2vw -i $fastaBatch -t $taxidBatch -k $K -d $dico | awk 'BEGIN{srand($SEED);} {printf "%06d %s\n", rand()*1000000, $0;}' | sort -n -T $TMPDIR | cut -c8- | vw --random_seed $SEED -f ${modelPrefix}_batch-${i}.model.sr -i ${modelPrefix}_batch-$(($i - 1)).model.sr --save_resume --progress 100000

		# last iteration : build final model (no save-resume mechanism)
		if [[ $i -eq ${NBATCHES} ]]; then
			$fasta2vw -i $fastaBatch -t $taxidBatch -k $K -d $dico | awk 'BEGIN{srand($SEED);} {printf "%06d %s\n", rand()*1000000, $0;}' | sort -n -T $TMPDIR | cut -c8- | vw --random_seed $SEED -f ${modelPrefix}_batch-${i}.model -i ${modelPrefix}_batch-$(($i - 1)).model.sr --progress 100000
		fi

		# remove old model
		rm ${modelPrefix}_batch-$(($i - 1)).model.sr
	fi

	# remove temp files
	rm  $outputDir/train.batch-$i*
done



