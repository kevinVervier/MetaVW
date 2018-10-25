

# specify (binary) model
model=../../2-build-models/output/train_small-db/vw-model_batch-10.bin

# specify test data
L=200
fasta=../../../data/small/test-dataset/fragments_L-$L.fasta

# specify output file
preds=../output/preds.fragments_L-$L.from-spectrum

# make prediction
spectrumpredict $model $fasta $preds

