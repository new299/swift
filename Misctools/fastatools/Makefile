RM = /bin/rm -f

C++ = g++
CCFLAGS = -O3 -g -I. -I../fastqtools

LDFLAGS = 

.PHONY: clean

all: fasta2flat fasta2fastq

fasta2flat: fasta2flat.cpp fastareader.h sequence.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) fasta2flat.cpp -o fasta2flat $(LDFLAGS)

fasta2fastq: fasta2fastq.cpp fastareader.h sequence.h ../fastqtools/fastqwriter.h ../fastqtools/scoredsequence.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) fasta2fastq.cpp -o fasta2fastq $(LDFLAGS) 

install: fasta2flat
	@echo -- $@ ---
	cp fasta2flat ../bin

clean:
	@echo --- $@ ---
	-$(RM) *.o core
