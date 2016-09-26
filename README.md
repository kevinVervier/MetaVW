# MetaVW: large-scale machine learning for metagenomics sequence classification
Source for reproducing results presented in [Large-scale machine learning for metagenomics sequence classification (Vervier et al., 2016)] (http://bioinformatics.oxfordjournals.org/content/32/7/1023.full.pdf+html)

## Data and code used in this project

The following file contains the data used in the paper (the small and large databases), as well as programs allowing to reproduce the results:
* [large-scale-metagenomics-1.0.tar.gz] (http://cbio.mines-paristech.fr/largescalemetagenomics/large-scale-metagenomics-1.0.tar.gz) (2.5 Go)

This archive is structured as follows :

* the _data_ directory contains the sequence data used to carry out the experiments
* the _src_ directory contains script allowing to reproduce the experiments
* the _tools_ directory contains source code of two utilities used to draw fragments from genomic sequences

We now detail the contents of the three directories present in the archive.
### data directory

* test-dataset/
  * test-db.abundance: a two-columns file (sequence_id : relative_abundance) required for simulations with Grinder software.
  * test-db.fasta: a FASTA file that contains the 398 sequences related to 193 species used to simulate validation sets.
  * test-db.meta: a five-columns file containing metadata on test-db genomes, like sequence.id, genome.id and the corresponding taxon IDs at strain, species and genus levels.
  * test-db.species-level.taxid: a one-column file that contains the species-level taxid for each sequence in the _test_ database.

* train-dataset/
  * train_large-db.fasta: a FASTA file that contains the 2768 sequences related to 774 species used to train on the _large_ database.
  * train_large-db.meta: a five-columns file containing metadata on large-db genomes, like sequence.id, genome.id and the corresponding taxon IDs at strain, species and genus levels.
  * train_large-db.species-level.taxid: a one-column file that contains the species-level taxid for each sequence in the _large_ database.
  * train_small-db.fasta: a FASTA file that contains the 1564 sequences related to 193 species used to train on the _small_ database.
  * train_small-db.meta: a five-columns file containing metadata on small-db genomes, like sequence.id, genome.id and the corresponding taxon IDs at strain, species and genus levels.
  * train_small-db.species-level.taxid: a one-column file that contains the species-level taxid for each sequence in the _small_ database.

### src directory

Each _src_ directory contains three sub-directories (input, output and src) allowing to reproduce the results presented in our paper (see the Tutorial section below):

* 1-generate-test-datasets/
  * input/
    * grinder.profile: contains parameters for Grinder software.
    * grinder-b-parameter_vs_mut-freq: a two-columns file providing equivalence between Grinder parameter for mutation rate and observed mutation rate.
  * output/
     * All the simulated data (fragments and sequences) will be stored in this directory.
  * src/
    * 01.generate-dataset-fragments.sh: Use drawfrag tool to generate fragments of length L from the test database.
    * 02.generate-dataset-reads-homo.sh: Use Grinder software to simulate reads, including homopolymer errors, from the test database.
    * 03.generate-dataset-reads-mutation.sh: Use Grinder software to simulate reads, including mutation errors, from the test database.
    * convert-gi-to-taxid.R: Match genome.id on species-level taxid.
* 2-build-models/
  * output/
    * All the sequences converted into VW bag-of-words will be (temporarily) stored in this directory, and discarded after learning.
    * The resulting classification models will also be stored in this directory. 
  * src/
    * 01.main.sh: Use VW to iteratively learn models.
* 3-make-predictions/
  * output/
    * All the VW predictions will be stored in this directory.
  * src/
    * 01.make-predictions.sh: Use VW to predict labels on validation sets.
    * 02.generate-graphs.R: Create result graphics similar to those present in the article.
    * vw-class-to-taxid.R: Match VW class labels on species-level taxid.

### tools directory

This directory contains source codes of two C programs involved in the Vowpal Wabbit (VW) based read classification pipeline :

* drawfrag is used to draw fragments from reference genomes
* fasta2vw is used to convert a fasta file into a plain text file compliant with VW

## Third-party softwares

To run the experiments, you need in addition to install the following third-party softwares:

* [Vowpal Wabbit] (https://github.com/JohnLangford/vowpal_wabbit/wiki) (version 7.7.0), the machine learning algorithm
* [Grinder] (http://sourceforge.net/projects/biogrinder/) (version 0.5.3) to simulate datasets
* [R] (http://cran.r-project.org/) (version 3.0.0) to analyze and plot the results

Optionally, you may want to install also the following softwares which we used in our paper for comparison to existing methods (but you do not need them if you just want to run our method):

* [FCP](http://kiwi.cs.dal.ca/Software/FCP) (version 1.0.6), the Fragment Classification Package which provides an implementation of Naive Bayes
* [BWA-MEM] (http://bio-bwa.sourceforge.net/) (version 0.7.4), the alignment-based approach tested in our benchmark

Our implementation uses the following libraries, which you do not need to install since they are provided in our software for convenience

* [kseq] (http://lh3lh3.users.sourceforge.net/kseq.shtml) library to parse FASTA/FASTQ files
* [GDL] (http://eqtnminer.sourceforge.net/index.php/Genetic_Data_analysis_Library_(GDL)) library implementing several structures like lists and hashtables (note that we provide a slightly updated version that properly compiles on Mac OS X)

# Tutorial: reproducing the results of the paper

After downloading the [project archive] (http://cbio.mines-paristech.fr/largescalemetagenomics/large-scale-metagenomics-1.0.tar.gz), first untar it:

    tar zxvf large-scale-metagenomics-1.0.tar.gz

Then, run the BASH script INSTALL.sh that can be found in the _tools_ directory. This script will process the installation of the GDL library, and create the binary executables:

    cd large-scale-metagenomics-1.0/tools
    sh INSTALL.sh

In order to check that everything went well during the installation, please use test.sh in the _tools/test_ directory. It will use installed tools to simulate a small dataset.

    cd test/
    sh test.sh
    ls output/
You should see the following files:    
_frags.fasta	frags.gi2taxid	frags.taxid	frags.vw	vw-dico.txt_

The following instructions have to be executed in the given order. At each step, we also point on tunable parameters.

1. ```cd src/1-generate-test-datasets/src```

2. ```sh 01.generate-dataset-fragments.sh```
    Generates the fragments test dataset based on the following parameters:
  * L: fragment length (default: 200)
  * COVERAGE: mean coverage for all genomes in the database (default: 1)

3. ``` sh 02.generate-dataset-reads-homo.sh```
  Generates the homopolymer test datasets based on the following parameters:
  * L: fragment length (default: 200)
  * COVERAGE: mean coverage for all genomes in the database (default: 1)
  * model: homopolymer errors (Balzer, Margulies, Richter)

4. ``` sh 03.generate-dataset-reads-mutation.sh```
  Generates the mutations test datasets based on the following parameters:
  * L: fragment length (default: 200)
  * COVERAGE: mean coverage for all genomes in the database (default: 1)

5. ``` cd ../../2-build-models/src```
6. ```sh 01.main.sh```
  Uses VW for an iterative learning (WARNING: can take time!) based on the following parameters:
 * DB: sequences database (small or large)
 * NBATCHES: number of train batches (default: 10)
 * L: fragment length (default: 200)
 * COVERAGE: mean coverage for all genomes in the database
 * K: k-mer size (default: 12)
 * LAMBDA1: L1 regularization parameter (default: 0)
 * LAMBDA2: L2 regularization parameter (default: 0)
 * NPASSES: number of training passes on each batch (default: 1)
 * BITS: hash table size (default: 31)

*WARNING* : this process takes some time, generates large files (~ 12 and 25 GB for the small and large databases, respectively) and has a comparable memory footprint.

7. ``` cd ../../3-make-predictions/src```

8. ```$ sh 01.make-predictions.sh```
Uses VW for predictions on validation sets based on the following parameters :
   * DB: sequences database (small or large)
   * NBATCHES: number of train batches (default: 10)
   * L: fragment length (default: 200)
   * K: k-mer size (default: 12)

9. ``` Rscript 02.generate-graphs.R```
Uses the previous results to plot the performance indicator considered (median species-level accuray), based on the following parameters:
* DB: sequences database (small or large)
* NBATCHES: number of train batches (default: 10)

# Classifying your own metagenomics sequences

You can train your own classifier by using the same scripts but changing some of the parameters.

In particular, you should be able to easily adapt the script 2-build-models/src/01.main.sh

Please remember that some parameters are mandatory:

* a training database (genome sequences stored as a multi-fasta file) and the corresponding taxids (defining the classes of the classification problem)
* a number of training batches to consider (NBATCHES)
* a mean coverage value to considere while generating training batches (COVERAGE)
* a k-mer size (K)

You may also want to directly use our models to make predictions for your own sequencing data (require training first).

To do this, please change the scripts in 3-make-predictions/src/ in order to point to your FASTA file.
