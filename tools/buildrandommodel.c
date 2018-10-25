
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#define ALPHABETSIZE 4  // DNA alphabet size


int power(int base, int exp)
// Simple power function for integer
{
    int result = 1;
    while(exp) { result *= base; exp--; }
    return result;
}



int main (int argc, char *argv[])
{    
	size_t nclasses, k, length, i, j;
	float w;
	FILE *fout;
  
    	// read parameters
	if (argc != 4) {
		fprintf(stderr, "Usage: %s nclasses k outputfile\n", argv[0]);
		return 1;
	}
	nclasses = atoi(argv[1]);
	k = atoi(argv[2]);
	fout = fopen(argv[3], "w");

	// compute number of features
	length = power(ALPHABETSIZE, k);

	// write number of classes and kmer length
	fwrite(&nclasses, sizeof(size_t), 1, fout);
	fwrite(&k, sizeof(size_t), 1, fout);

	// write model values
	for(i = 0; i < length; i++){		// NB : respect format of spectrumpredict (feature-level 'blocs' of nclasses variables)
		for(j = 0; j < nclasses; j++){
			w = (float)(rand())/(float)(RAND_MAX);
			fwrite(&w, sizeof(float), 1, fout);
		}
	}

	// write intercepts
	for(i = 0; i < nclasses; i++){
		w = (float)(rand())/(float)(RAND_MAX);
		fwrite(&w, sizeof(float), 1, fout);
	}

	// close output file
	fclose(fout);

	return 0;
}


