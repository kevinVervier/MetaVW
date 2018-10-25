

GDLDIR=./ext/gdl-1.2.1/GDL/

CC=gcc
CFLAGS= -O2 -Wall -Wshadow -I./ -I${GDLDIR}/include -c
OBJS = drawfrag.o fasta2vw.o enumerateKmers.o parseHash.o buildrandommodel.o  spectrumpredict.o  spectrumpredict_scores.o

all: drawfrag fasta2vw enumerateKmers parseHash buildrandommodel spectrumpredict spectrumpredict_scores

.c.o :
	$(CC) $(CFLAGS) $<

drawfrag : drawfrag.o
	$(CC) -o $@ drawfrag.o -lz -L${GDLDIR}/lib -lgdl -lm

fasta2vw : fasta2vw.o
	$(CC) -o $@ fasta2vw.o -lz -L${GDLDIR}/lib -lgdl -lm

enumerateKmers : enumerateKmers.o
	$(CC) -o $@ enumerateKmers.o -lz -L${GDLDIR}/lib -lgdl -lm

parseHash : parseHash.o
	$(CC) -o $@ parseHash.o -lz -L${GDLDIR}/lib -lgdl -lm

buildrandommodel : buildrandommodel.o
	$(CC) -o $@  buildrandommodel.o

spectrumpredict : spectrumpredict.o
	$(CC) -o $@ spectrumpredict.o -lz

spectrumpredict_scores : spectrumpredict_scores.o
	$(CC) -o $@ spectrumpredict_scores.o -lz


clean :
	rm $(OBJS) drawfrag fasta2vw enumerateKmers parseHash buildrandommodel spectrumpredict spectrumpredict_scores

