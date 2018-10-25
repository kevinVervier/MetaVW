#!/bin/sh

# read paramers
K=$1
N=$2
model=$3
hash=$4
hashBin=$5

echo "converting vw model to binary"
# enumerate kmers
echo "    -1 - enumerating kmers of length $K"
allKmers=TMP_all-kmers_k-$K.vw
enumerateKmers -k $K -o $allKmers

# invert hash table
echo "    -2 - inverting vw hash table"
vw -i $model -d $allKmers --invert_hash $hash

# convert hash to binary
echo "    -3 - converting hash to binary"
parseHash -i $hash -o $hashBin -k $K -n $N

# remove temp file
#rm $allKmers
#rm $hash


