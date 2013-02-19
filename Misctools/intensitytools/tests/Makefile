RM = /bin/rm -f

C++ = g++
CCFLAGS = -I. -I..

LDFLAGS = 

.PHONY: test

all: test

test: testmain.cpp utf.h ../*.h ../*.cpp ./*.cpp ./*.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) testmain.cpp test_intensity2filterblobs.cpp test_intensity2pairwiseplot.cpp test_readintensity.cpp test_intensitysequence.cpp -o test $(LDFLAGS)

clean:
	@echo --- $@ ---
	-$(RM) *.o core
