
# specify K
K=10


# get 1st sequence
L=200
fasta=../../../data/small/test-dataset/fragments_L-$L.fasta
head -n 4 $fasta > DEBUG_1st-seq.fasta

# extract kmers
fasta2vw -i  DEBUG_1st-seq.fasta -k $K > DEBUG_1st-seq.fasta2vw

# extract weights of 11th class
hash=../../2-build-models/output/train_small-db/vw-model_batch-10.hash
grep "\[10\]" $hash > DEBUG_weights

# compute score
Rscript 04.DEBUG_computeScore.R DEBUG_1st-seq.fasta2vw DEBUG_weights

# check with vw score
head -n 1 ../output/preds.fragments_L-200.from-vw.scores
head -n 1 ../output/preds.fragments_L-200.from-vw.scores | cut -f 11 -d ' '

# check with score compute by spectrumpredict
model=../../2-build-models/output/train_small-db/vw-model_batch-10.bin
spectrumpredict_scores $model DEBUG_1st-seq.fasta tt


