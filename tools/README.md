This directory contains source codes of two C programs involved in the Vowpam Wabbit (VW) based read classification pipeline :
* `drawfrag' is used to draw fragments from reference genomes
* `fasta2vw' is used to convert a fasta file into a plain text file compliant with VW

These two tools rely on two existing libraries :
* the GDL library that implements (among  other thing) many file manipulation utilities and data structures (e.g., list and hash tables)
* the KSEQ library that is used to manipulate fasta and fastq files
These two libraries are included in the "ext" directory.

The `INSTALL.sh' bash script carries out the installation of these tools : 
* installation of the GDL library (installed here in ext/gdl-1.2/GDL)
* installation of drawfrag and fasta2vw.

To run is, you should type at the prompt:

``` sh INSTALL.sh ```

The `test' directory contains a simple script illustrating how to use these tools on a small fasta file.
