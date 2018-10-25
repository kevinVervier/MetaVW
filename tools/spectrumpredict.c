//
//  spectrumpredict.c
//  
//
//  Created by Jean-Philippe Vert on 04/12/2014.
//
//  Fast scoring of a DNA sequence with a model that adds the weights of all the k-mers contained in the sequence

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <zlib.h>
#include <ext/kseq.h>
KSEQ_INIT(gzFile, gzread);


#define MAXK 50     // Max k-mer use
#define ALPHABETSIZE 4  // DNA alphabet size

typedef struct {
    size_t nclasses, k, length;
    float *weight;
} weight_t;
// A structure to store the vector of weights. k is the kmer size; length the length of the weight vector (should be power(ALPHABETSIZE,k). weight is the vector of weights


int power(int base, int exp)
// Simple power function for integer
{
    int result = 1;
    while(exp) { result *= base; exp--; }
    return result;
}

unsigned int DNAtoint(char c)
// Convert a nucleotide char to an integer from 0 to 3
{
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
//            fprintf(stderr,"Error: %c is not a nucleotide\n",c);
//            exit(2);
            return -1;
    }
}

unsigned int dnaseqtoint(char *dnaseq)
// Transform a DNA sequence to an int (by concatenating the 2-bit representations of the 4 bases)
{
    unsigned int x=0,i;
    for (i=0; i<strlen(dnaseq); x=(x<<2)|DNAtoint(dnaseq[i++]));
    return x;
}


weight_t *readmodelfromfile(char *filename)
// Read a weight vector from a file
// Each row of the file should contain a kmer of DNA and a weight, separated by a TAB
// The k-mer length is inferred from the first k-mer and the length of the weight vector set accordingly
{
    FILE *inp;
    int first=1;
    char s[MAXK];
    double d;
    weight_t *w = malloc(sizeof(weight_t));


    /* Open file */
    if ((inp=fopen(filename,"r")) == NULL) {
        fprintf(stderr,"Error: cannot open file %s\n",filename);
        exit(2);
    }
    
    while (fscanf(inp,"%s\t%lf",s,&d) ==2) {
        if (first) {
            // Estimate k-mer size and allocate memory for the weight vector
            w->k = strlen(s);
            w->length = power(ALPHABETSIZE,w->k);
            if ((w->weight=(float *)malloc(sizeof(float)*w->length)) == NULL) {
                fprintf(stderr,"Error: run out of memory\n");
                exit(2);
            }
            first=0;
        }
        // Store the weight
        w->weight[dnaseqtoint(s)] = d;

    }
    
    fclose(inp);
    return w;
}


weight_t *readmodelfrombinaryfile_fast(char *filename)
// Read a weight vector from a binary file
// the 1st value contains the number of classes
// the 2nd value contains the kmer size
// the  [nclasses x 4^(k+1)] following values contain the model parameters 
{
	FILE *inp;
	int ltotal;
	weight_t *w = malloc(sizeof(weight_t));

	// Open file
	if ((inp=fopen(filename,"r")) == NULL) {
        fprintf(stderr,"Error: cannot open file %s\n",filename);
        exit(2);
	}
    
	// Read number of classes and kmer length
	fread(&(w->nclasses), sizeof(size_t), 1, inp);
	fread(&(w->k), sizeof(size_t), 1, inp);

	// Compute number of features
	w->length = power(ALPHABETSIZE,w->k);
    
	// allocate memory for weight vector
	ltotal = w->nclasses * (w->length + 1); // NB : add 1 to feature counts to store intercepts
	w->weight = (float *)malloc(sizeof(float)*ltotal);
	
	// read weights
	fread(w->weight, sizeof(float), ltotal, inp);

	// Close input file
	fclose(inp);
	// Return
	return w;
}



/*

static int kseq_computescore(kseq_t *seq, weight_t *w, float *score)									\
{																	\
    int c,k,j;
    unsigned int x=0,mask=w->length-1;
    kstream_t *ks = seq->f;											\
    size_t topclass;
    float scoremax;


    for(k=0;k<w->nclasses;k++)
        score[k] = 0;

    if (seq->last_char == 0) {  \
        while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@');	\
        if (c == -1) return -1; 				\
        seq->last_char = c;											\
    } 				\
//    seq->comment.l = seq->seq.l = seq->qual.l = 0;
//    if (ks_getuntil(ks, 0, &seq->name, &c) < 0) return -1;			\
    
    // Skip the header line
    while ((c = ks_getc(ks)) != -1 && c != '\n');

    // Read the sequence
    
    // Read the first k-mer
    for (j=0; j<w->k ;) {
//        x=(x<<2)|DNAtoint(ks_getc(ks));
        c=DNAtoint(ks_getc(ks));
        if (c>=0) {
            x=(x<<2)|c;
            j++;
        }
    }
    
    // Store the scores of the first kmer
    for(k=0; k < w->nclasses; k++)
            score[k] += w->weight[k][x];
    // Read subsequent letters
    while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') {
        j=DNAtoint(c);
        if (j>=0) {
            x=((x<<2)& mask)|j;
            for(k=0; k < w->nclasses; k++)
                score[k] += w->weight[k][x];
        }
    }
    if (c == '>' || c == '@') seq->last_char = c; // the first header char has been read 	\
    
        
        
    // Compute best class
    topclass = 0;
    scoremax = score[0];
    for(k = 1; k < w->nclasses; k++){
        if(score[k] > scoremax){
            scoremax = score[k];
            topclass = k;
        }
    }
    return topclass;
}
*/

static int kseq_computescore_fast(kseq_t *seq, weight_t *w, float *score)
{
    int c, k, j, nc=w->nclasses, nfeats=w->length, K = w->k;
    unsigned int x=0,mask=w->length-1;
    kstream_t *ks = seq->f;
    size_t topclass;
    float scoremax, *myw;   

    // set scores to 0
    for(k=0; k < nc; score[k++] = 0);  	

    // Skip the header line
    while ((c = ks_getc(ks)) != -1 && c != '\n');

// DEBUG : checking header
//    printf("\n\n");
//    printf("-> skipping header\n");
//    printf("\t- header =\n");
//    while ((c = ks_getc(ks)) != -1 && c != '\n'){
//	printf("%c", c);	    
//    };

    if (c==-1) return -1; //EOF


//    // Read the first k-mer
//    for (j=0; j<K ;) {
//  	c=DNAtoint(ks_getc(ks));
//        if (c>=0) {
//            x=(x<<2)|c;
//            j++;
//        }
//    }

    // Read the 1st kmer --> MODIFIED TO HANDLE CASES WHERE SEQUENCE MADE OF ONLY N's !
    int cc;
    for (j=0; j<K ;) {
	cc = ks_getc(ks);
	if(cc == '>'){	// WARNING: if we reach the next line, means no kmer found --> return class 0 (no prediction)
		return(0);
	}
  	c=DNAtoint(cc);
        if (c>=0) { 
            x=(x<<2)|c;
            j++;
        }
    }

    // Store the scores of the first kmer
    myw=w->weight + x*nc;
    for(k=0; k < nc; score[k++] = *myw++);

    // Read subsequent letters
    while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') {
        j=DNAtoint(c);
        if (j>=0) {
            x=((x<<2)& mask)|j;
            myw=w->weight + x*nc;	// NB : weights are organized as feature-level blocs (1st set of nclass values = feature 1 and so on) 
            for(k=0; k < nc; k++) score[k] += 2*myw[k];	// NB : multiply by 2 for RT kmer (NB : checked, weights are similar)
        }
    }

    // add intercepts
    myw=w->weight + nfeats*nc;		// NB : intercepts are stored at the end of the vector (as an additional set of nclass values)
    for(k=0; k < nc; k++) score[k] += myw[k];

    // compute best class
    topclass = 0;
    scoremax = score[0];
    for(k = 1; k < nc; k++){
        if(score[k] > scoremax){
            scoremax = score[k];
            topclass = k;
        }
    }
    topclass = topclass+1;


//    // DEBUG : print scores
//    for(k = 0; k < nc; k++){
//	printf("%.6f ", score[k]);
//    }
//    printf("\n");

    return topclass;
}




int main (int argc, char *argv[])
{
    
	weight_t *w;
	gzFile fin;
	FILE *fout;
	kseq_t  *seq;
	time_t start, end;
	double diff;
	int i;

	if (argc != 4) {
		fprintf(stderr, "Usage: %s model fastafile outputfile\n", argv[0]);
		return 1;
	}
    
	// Read the vector of weights (of each kmer)
	// Note that if a kmer is not in the file, it is just assigned a weight of 0
	time(&start);
	w = readmodelfrombinaryfile_fast(argv[1]);
	time(&end);
	diff = difftime (end,start);
	printf("reading model involving %d classes with %d features (k = %d) took %.2lf seconds\n", w->nclasses, w->length, w->k, diff);

 	// open input file
	fin = gzopen(argv[2], "r");
	seq = kseq_init(fin);

	//  open output file
	fout = fopen(argv[3], "w");

	// Allocate memory once for score
	float *score = (float *)malloc(sizeof(float)*w->nclasses);

	// process sequences
	int cptr = 0;
	time(&start);	

	//while ((i=kseq_computescore(seq,w,score)) >= 0){
	while ((i=kseq_computescore_fast(seq,w, score)) >= 0){
		cptr++;
        	if(cptr > 0 && cptr % 100000 == 0){
			time(&end);
			diff = difftime (end,start);
			printf("\t- processing sequence no %d\n", cptr);
			printf("\t\t--> last loop took %.2f seconds\n", diff);
			time(&start);
		}			
		fprintf(fout, "%lu\n",i);
	}

	// close files
	fclose(fout);
	gzclose(fin);
 
     	// free memory
	free(score);
	free(w->weight);
	free(w);
	kseq_destroy(seq);

	// exit
	return 1;
}



