#Copyright (c) by respective owners including bioMerieux, and
#individual contributors Kevin Vervier and Pierre Mahe. All #rights reserved.  Released under a BSD (revised)
#license as described in the file LICENSE.


# specify input  data
fasta=../../../data/medium-large/test-dataset/test_db.fasta
taxids=../../../data/medium-large/test-dataset/test_db.species-level.taxid
abundance=../../../data/medium-large/test-dataset/test_db.abundance
meta=../../../data/medium-large/test-dataset/test_db.meta

# specify fragments parameters (length and cover ~ number)
size=200
cover=1


#############################
###### SMALL SCALE TEST #####
cover=0.05
#############################


# set seed (for reproducibility)
seed=42

# draw reads
for model in Balzer Margulies Richter
do
	echo generating reads with model $model 
	# call grinder
	grinder -bn test.homo-$model -cf $cover -rd $size -rs $seed -pf ../input/grinder.profile -af $abundance -od ../output -rf $fasta -hd $model
	# extract genome ids from fastq headers
	grep '@' ../output/test.homo-$model-reads.fastq | cut -f 2 -d ' ' | cut -f 2 -d '=' > ../output/test.homo-$model-reads.gis
	# convert to taxids
	Rscript convert-gi-to-taxid.R ../output/test.homo-$model-reads.gis $meta ../output/test.homo-$model-reads.taxid
done


