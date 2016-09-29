

GDLDIR=./ext/gdl-1.2/GDL
#GDLDIR=/Users/jp/tmp/gdl/GDL/

CC=gcc
CFLAGS= -O2 -Wall -Wshadow -Wextra -I./ -I${GDLDIR}/include -c
OBJS = drawfrag.o fasta2vw.o

all: drawfrag fasta2vw

.c.o :
	$(CC) $(CFLAGS) $<

drawfrag : drawfrag.o
	$(CC) -o $@ drawfrag.o -lz -L${GDLDIR}/lib -lgdl -lgsl -lm

fasta2vw : fasta2vw.o
	$(CC) -o $@ fasta2vw.o -lz -L${GDLDIR}/lib -lgdl -lgsl -lm

clean :
	rm $(OBJS) drawfrag fasta2vw
