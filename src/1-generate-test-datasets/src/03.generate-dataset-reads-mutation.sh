#Copyright (c) by respective owners including bioMerieux, and
#individual contributors Kevin Vervier and Pierre Mahe. All #rights reserved.  Released under a BSD (revised)
#license as described in the file LICENSE.

# specify input  data
fasta=../../../data/test-dataset/test_db.fasta
taxids=../../../data/test-dataset/test_db.species-level.taxid
abundance=../../../data/test-dataset/test_db.abundance
meta=../../../data/test-dataset/test_db.meta

# specify fragments parameters (length and cover ~ number)
size=200
cover=1

#############################
###### SMALL SCALE TEST #####
#cover=0.05
#############################


# set seed (for reproducibility)
seed=42

# draw reads
while read myLine
do
	# extract b parameter and corresponding mutation frequency
	coeff=$(echo $myLine | cut -f 1 -d ' ')
	freq=$(echo $myLine | cut -f 2 -d ' ')

	# call grinder
	grinder -bn test.mutation-$freq -cf $cover -rd $size -rs $seed -pf ../input/grinder.profile -af $abundance -od ../output -rf $fasta -md poly4 3e-3 $coeff
	# extract genome ids from fastq headers
	grep '@' ../output/test.mutation-$freq-reads.fastq | cut -f 2 -d ' ' | cut -f 2 -d '=' > ../output/test.mutation-$freq-reads.gis
	# convert to taxids
	Rscript convert-gi-to-taxid.R ../output/test.mutation-$freq-reads.gis $meta ../output/test.mutation-$freq-reads.taxid

done < ../input/grinder-b-parameter_vs_mut-freq.txt



