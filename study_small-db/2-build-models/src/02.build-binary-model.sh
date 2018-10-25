
# specify parameters
K=10
N=51
model=../output/train_small-db/vw-model_batch-10.model
hash=../output/train_small-db/vw-model_batch-10.hash	# NB : this file can get pretty big (but is removed in the end)
hashBin=../output/train_small-db/vw-model_batch-10.bin

vw-to-binary.sh $K $N $model $hash $hashBin

