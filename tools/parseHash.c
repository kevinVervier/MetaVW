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

static gdl_string * PROGRAM = "parseHash";

static int help_flag    = 0;
static int verbose_flag = 0;
static int debug_flag = 0;

static gdl_string * INPUT = NULL;
static gdl_string * OUTPUT = NULL;
static int K = 0;
static int N = 0;

static struct option long_options[] =
{
		/* These options set a flag. */
		{"help", no_argument,		&help_flag, 1},
		{"verbose", no_argument,	&verbose_flag, 1},
		{"debug", no_argument,		&debug_flag, 1},
		/* These options don't set a flag.
         	We distinguish them by their indices. */
		{"input",   required_argument, 0, 'i'},
		{"output",   required_argument, 0, 'o'},
		{"kmer", required_argument, 0, 'k'},
		{"nclasses", required_argument, 0, 'n'},
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

		c = getopt_long (argc, argv, "i:o:k:n:",
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
		case 'i':
			INPUT = gdl_string_clone (optarg);
			break;
		case 'o':
			OUTPUT = gdl_string_clone (optarg);
			break;
		case 'k':
			K = atoi(optarg);
			break;
		case 'n':
			N = atoi(optarg);
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
	if(N == 0){
		GDL_ERROR_VAL ("No number of classes provided", GDL_FAILURE, 1);
	}
	if(INPUT == 0){
		GDL_ERROR_VAL ("No input file provided", GDL_FAILURE, 1);
	}
	if(OUTPUT == 0){
		GDL_ERROR_VAL ("No output file provided", GDL_FAILURE, 1);
	}
	return GDL_SUCCESS;
}

static int
help (void)
{
	printf("%s - version 1.0\n", PROGRAM);
   	printf("Copyright (C) 2017 bioMerieux, France\n");
	printf ("\n");
	printf("Transforms a (hash-inverted) Vowpal Wabbit (VW) model into a binary file\n");
	printf ("\n");
	printf ("--help\tDisplay a brief help on program usage\n");
	printf ("--verbose\tOutput message on standard output to see what the program is doing\n");
	printf ("--debug\t debug mode; output many messages on standard output to see what the program is doing\n");
	printf ("\n");
	printf ("--input or -i\t VW model in plain text format, obtained using invert_hash option\n");
	printf ("--output or -o\t output file\n");
	printf ("--kmer or -k\t k-mers length\n");
	printf ("--nclasses or -n\t number of classes in the model\n");
	return GDL_SUCCESS;
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
	
		FILE *stream_in = 0, *stream_out = 0;
		size_t nkmers, weight_length, n, ntok, ntok2, class, x, j, mask, index;
		int m;
		float *weight, w;
		gdl_string    * line = 0;
		gdl_string ** toks, **toks2;

		// set verbose to true if debug
		if(debug_flag){
			verbose_flag = 1;
		}

		// initialize vector of weights
		nkmers = power(ALPHABETSIZE,K);
		mask = nkmers-1; // use to convert kmers to integers
		weight_length = N * (nkmers+1);
	       	weight = (float *)calloc(weight_length, sizeof(float));

		// process data 
		stream_in = gdl_fileopen (INPUT, "r");
		while(gdl_getline (&line, &n, stream_in)!=-1){
			// split according to :
			ntok = 0;
			toks = gdl_string_split(line, ":", &ntok);
			if(ntok == 3){ // if fewer than 3 elements : header
				// extract weight
				w = atof(toks[2]);

				// extract kmer & class
					// split 1st field according to [
				ntok2 = 0;
				toks2 = gdl_string_split(toks[0], "[", &ntok2);
					// if no [] found : corresponds to 1st class
				if(ntok2 == 1){
					class = 1;
					if(debug_flag){
						printf("%s : kmer %s of class no %d with weight %f\n", line, toks2[0], class, w);
					}
				}
					// else : extract class
				if(ntok2 == 2){
					// discard lat character
					toks2[1][strlen(toks2[1])-1] = '\0';
					class = atoi(toks2[1]) + 1;
					if(debug_flag){
						printf("%s : kmer %s of class no %d with weight %f\n", line, toks2[0], class, w);
					}
				}
				// --> at this point we have : weight (w), class id (class) and kmer sequence (toks2[0]

				// get index corresponding to the kmer
				if(strcmp(toks2[0], "Constant") == 0){ // handling of intercept : last feature
					x = nkmers;
					if(debug_flag){
						printf("%s\n", line);
					}
				}else{
					x = 0;
					for(j = 0; j < strlen(toks2[0]); j++){
						toks2[0][j] = toupper(toks2[0][j]);
						m = DNAtoint(toks2[0][j]);
						if (m>=0) {
							 x = ((x<<2)& mask)|m;
						}else{
							 x = -1;
						}
					}
				}
				// store value in vector or weights
				if(debug_flag){
					printf("%s : kmer %s (index %d) of class no %d with weight %f\n", line, toks2[0], x, class, w);
				}
				index = (class - 1) + x*N;
				weight[index] = w;
			}

		    	gdl_string_free (line);
		    	line=0;
		}
    	// write output file
	stream_out = gdl_fileopen (OUTPUT, "w");
	fwrite(&N, sizeof(size_t), 1, stream_out);
	fwrite(&K, sizeof(size_t), 1, stream_out);
	fwrite(weight, sizeof(float), weight_length, stream_out);

	// close files
	gdl_fileclose(INPUT, stream_in);
	gdl_fileclose(OUTPUT, stream_out);
	}
	exit (0);
}


