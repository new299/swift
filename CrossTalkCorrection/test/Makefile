
CPPFILES = 

all:
	g++ testmain.cpp test_crosstalkcorrection.cpp $(CPPFILES) -I.. -I../../include -I../../lib/gslmm/include -I../../lib/gsl/include `../../lib/gsl/bin/gsl-config --libs` -g -DHAVE_FFTW $(CPPFLAGS) -O3 -fomit-frame-pointer -DFTYPE=float -Wall -Wsign-compare -Wpointer-arith -o test
