# draw fragments of size 20 covering each sequence 0.8 times
../drawfrag -i input/seq.fasta -t input/seq.taxid -l 20 -c 0.8 -o output/frags.fasta -g output/frags.gi2taxid

# extract fragments taxids
cut -f 2 output/frags.gi2taxid > output/frags.taxid

# convert fragments to k-mers of size 6 in a format compliant with Vowpal Wabbit 
../fasta2vw -i output/frags.fasta -t output/frags.taxid -k 6 -o output/frags.vw -d output/vw-dico.txt



