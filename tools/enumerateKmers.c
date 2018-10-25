/*
Copyright (c) by respective owners including bioMerieux, and
individual contributors Kevin Vervier and Pierre Mahe. All rights reserved.  Released under a BSD (revised)
license as described in the file LICENSE.
 */
#include <getopt.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>
#include <math.h>

#include <gdl/gdl_common.h>
#include <gdl/gdl_version.h>
#include <gdl/gdl_errno.h>
#include <gdl/gdl_sort_uint.h>
#include <gdl/gdl_io.h>
#include <gdl/gdl_hash.h>
#include <gdl/gdl_list.h>
#include <gdl/gdl_string.h>
#include <gdl/gdl_util.h>

#include <zlib.h>
//#include <ext/kseq.h>
//KSEQ_INIT(gzFile, gzread);

#define ALPHABETSIZE 4  // DNA alphabet size

static gdl_string * PROGRAM = "enumerateKmers";

static int help_flag    = 0;
static int verbose_flag = 0;
static int stdout_flag = 1;
static int debug_flag = 0;

static gdl_string * OUTPUT = NULL;
//static gdl_string * K = NULL;
static int K = 0;

static struct option long_options[] =
{
		/* These options set a flag. */
		{"help", no_argument,		&help_flag, 1},
		{"verbose", no_argument,	&verbose_flag, 1},
		{"debug", no_argument,		&debug_flag, 1},
		/* These options don't set a flag.
         	We distinguish them by their indices. */
		{"output",   required_argument, 0, 'o'},
		{"kmer_size", required_argument, 0, 'k'},
		{0, 0, 0, 0}
};

static int
parse_argument (int argc, char *argv[])
{
	int c;

	while (1)
	{

		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "o:k:",
				long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
		case 0:
			/* If this option set a flag, do nothing else now. */
			if (long_options[option_index].flag != 0)
				break;
			printf ("option %s", long_options[option_index].name);
			if (optarg)
				printf (" with arg %s", optarg);
			printf ("\n");
			break;
		case 'o':
			OUTPUT = gdl_string_clone (optarg);
			stdout_flag = 0;
			break;
		case 'k':
			//K = gdl_string_clone (optarg);
			K = atoi(optarg);
			break;
		case '?':
			GDL_ERROR_VAL ("Unknown arguments", GDL_EINVAL, -1);
		default:
			GDL_ERROR_VAL ("Bad arguments", GDL_EINVAL, -1);
		}
	}
	return GDL_SUCCESS;
}

static int
check_argument (void)
{
	if(K == 0){
		GDL_ERROR_VAL ("No kmer size provided", GDL_FAILURE, 1);
	}
	return GDL_SUCCESS;
}

static int
help (void)
{
	printf("%s - version 1.0\n", PROGRAM);
   	printf("Copyright (C) 2017 bioMerieux, France\n");
	printf ("\n");
	printf("Enumerate (canonical) k-mers of a given length into a Vowpal Wabbit (VW) compliant output file\n");
	printf ("\n");
	printf ("--help\tDisplay a brief help on program usage\n");
	printf ("--verbose\tOutput message on standard output to see what the program is doing\n");
	printf ("--debug\t debug mode; output many messages on standard output to see what the program is doing\n");
	printf ("\n");
	printf ("--kmer or -k\t k-mer size\n");
	printf ("[--output or -o\t output file]\n");
	printf ("\t (if not specified, results printed on standard output)\n");
	return GDL_SUCCESS;
}


// compute reverse-complement of a k-mer
gdl_string * reverse_transcribe(gdl_string * s){
        // initialize result
        gdl_string * res;
        res = gdl_string_alloc(strlen(s));
        // process sequence
        size_t i;
        for (i = 0; i < strlen(s); i++){
                // reverse
                res[i] = s[strlen(s)-i-1];
                // transcribe
                switch(res[i]){
                        case 'A' : res[i] = 'T'; break;
                        case 'T' : res[i] = 'A'; break;
                        case 'G' : res[i] = 'C'; break;
                        case 'C' : res[i] = 'G'; break;
                        default : break;
                }
        }
        // return
        return(res);
}



// Simple power function for integer
int power(int base, int exp){
    int result = 1;
    while(exp) { result *= base; exp--; }
    return result;
}



// Convert a nucleotide char to an integer from 0 to 3
unsigned int DNAtoint(char c){
    switch (toupper(c)) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
		// fprintf(stderr,"Error: %c is not a nucleotide\n",c);
		// exit(2);
            return -1;
    }
}


// Convert an integer from 0 to 3 into a nucleotide
char intToDNA(int x){
    switch (x) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
		fprintf(stderr,"Error: %d is greater than 3\n",x);
		exit(2);
            return -1;
    }
}



int
main (int argc, char *argv[])
{

	int status;

	parse_argument (argc, argv);

	if (help_flag)
	{
		exit (help());
	}

	status = check_argument ();

	if (status == GDL_SUCCESS)
	{

	
		FILE * stream = 0;
		size_t nfeats, nfeats_cano, i, j, k, cptr, mask_char, mask_kmer, i_rt, cptr_written;
		int m;

		// set verbose to true if debug
		if(debug_flag){
			verbose_flag = 1;
		}

		// open output file
  		if(!stdout_flag)
		{
			stream = gdl_fileopen (OUTPUT, "w");
		}

		// define number of features
		nfeats = power(ALPHABETSIZE,K);		//4^K possible k-mers
		nfeats_cano = nfeats/2;                 // 4^K/2 canonical ones
		if(verbose_flag){
			printf("Enumerating kmers of length K = %d --> nfeats = %zd\n", K, nfeats); 
		}

		// define bits-masks
		mask_char = 3;
		mask_kmer = nfeats-1;                	


		// print VW heaer
		if(stdout_flag){
			printf("1 |");
		}else{
                 	fprintf(stream,"1 |");
		}

		// process each kmer
		cptr = 0;
		j = 0;
		cptr_written = 0;
		for(i = 0; i < nfeats; i++){
			if( (i > 0) & (i%100000 == 0) ){
				if(verbose_flag){
					printf("\n\n*** processing feature no %zd (out of %zd)***\n", i, nfeats);
				}
			}
			// init kmer
			gdl_string *kmer;
			kmer = 0;
			kmer = gdl_string_alloc (K);
			// build kmer from feature index
			cptr = i;
			j = cptr & mask_char;
			kmer[K-1] = intToDNA(j);
			for(k = 1; k < K ; k++){
				cptr = cptr >> 2;
				j = cptr & mask_char;
				kmer[K-k-1] = intToDNA(j);
			}
			if(debug_flag){
				printf("--> kmer sequence = %s (index = %zd) \n", kmer, i);
			}
			// print kmer
			cptr_written = cptr_written + 1;
			if(stdout_flag){
				printf(" %s", kmer);
			}else{
				fprintf(stream," %s", kmer);
			}
			// free kmers
			gdl_string_free(kmer);

			/* 
			// extension : print only canonical kmer - NB : works for odd valus of k only 
			// build RT-kmer
			gdl_string *rt_kmer;
			rt_kmer = 0;
			rt_kmer = reverse_transcribe(kmer);	
			// extract reverse index
			i_rt = 0;
			for(k = 0; k < K; k++){
				m = DNAtoint(rt_kmer[k]);
				i_rt = ((i_rt<<2)& mask_kmer)|m;
			}
			if(debug_flag){
				printf("--> rt-kmer sequence = %s (index = %zd) \n", rt_kmer, i_rt);
			}
			// only print kmer if canonical
			if(i < i_rt){
				cptr_written = cptr_written + 1;
				if(stdout_flag){
					printf(" %s", kmer);
				}else{
					fprintf(stream," %s", kmer);
				}
			}
			// free kmer
			gdl_string_free(rt_kmer);
			*/
		}

		// print end-of-line
		if(stdout_flag){
			printf("\n");
		}else{
			fprintf(stream,"\n");
		}

		// close output file
		if(!stdout_flag){
			gdl_fileclose (OUTPUT, stream);
		}

		// sanity check
		if(cptr_written != nfeats){
			printf("####### ERROR : not the good number of features (%zd instead of %zd expected) #######\n", cptr_written, nfeats); // NB : would not work if canonical only, need to compare against nfeats_cano
			exit(1);	
		}else{
			if(verbose_flag){
				printf("==> DONE : expected number of features written (%zd)\n", cptr_written);
			}
		}
	}
	exit (0);
}



