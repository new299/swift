RM = /bin/rm -f

C++ = g++
CCFLAGS = -O3 -g -I.

LDFLAGS = 

.PHONY: clean

all: fastq2duplicatecheck fastqfilter fastqascii2txtqual fastqtrim fastqpolyremove

fastqfilter: fastqfilter.cpp fastqreader.h fastqwriter.h scoredsequence.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) fastqfilter.cpp -o fastqfilter $(LDFLAGS)

fastqtrim: fastqtrim.cpp fastqreader.h fastqwriter.h scoredsequence.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) fastqtrim.cpp -o fastqtrim $(LDFLAGS)

fastqascii2txtqual: fastqascii2txtqual.cpp fastqreader.h fastqwriter.h scoredsequence.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) fastqascii2txtqual.cpp -o fastqascii2txt $(LDFLAGS)

fastq2split: fastq2split.cpp fastqreader.h fastqwriter.h scoredsequence.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) fastq2split.cpp -o fastq2split $(LDFLAGS)


fastqpolyremove: fastqployremove.cpp fastqreader.h fastqwriter.h scoredsequence.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) fastqployremove.cpp -o fastqpolyremove

fastq2duplicatecheck: fastq2duplicatecheck.cpp fastqreader.h scoredsequence.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) fastq2duplicatecheck.cpp -o fastq2duplicatecheck

install: fastqascii2txtqual fastqtrim fastqfilter
	cp fastqascii2txtqual ../bin
	cp fastqtrim ../bin
	cp fastqfilter ../bin

clean:
	@echo --- $@ ---
	-$(RM) *.o core


