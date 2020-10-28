SVNDEF := -D'SVN_REV="$(shell svnversion -n .)"'

SWIFT_IMAGEANALYSIS_H = ./SwiftImageAnalysis

SWIFT_IMAGEANALYSIS = $(SWIFT_IMAGEANALYSIS_H)/*.h \
                      $(SWIFT_IMAGEANALYSIS_H)/*.cpp \
                      ./PhasingCorrection/*.h ./PhasingCorrection/*.cpp ./CrossTalkCorrection/*.h ./CrossTalkCorrection/*.cpp ./Reporting/*.h ./SmallAlign/*.h

# TODO: Above list of dependencies is not complete

CPPFILES = ./include_lib/gnuplot_i.cc ./SwiftImageAnalysis/SwiftFFT.cpp 

HFILES = ./CrossTalkCorrection/*.cpp ./CrossTalkCorrection/*.h ./include/CommandLine.h

IFILES = -I./Filters -I./SmallAlign -I./MockImageAnalysis -I./BaseCaller \
         -I./include -I./CrossTalkCorrection -I./PhasingCorrection \
         -I./include_lib -I /software/solexa/include \
	 -I./Reporting

LFILES       = -L /software/solexa/lib -lgsl -lgslcblas -lfftw3f -ltiff
SWIFT_LFILES = -L /software/solexa/lib -lgsl -lgslcblas -lfftw3 -ltiff

CPPFLAGS = $(SVNDEF) -O3 -DHAVE_FFTW -DFTYPE=float -Wall -Wsign-compare -Wpointer-arith -std=c++14
#CPPFLAGS = -g -DHAVE_FFTW -Wpointer-arith
INTEL_CPPFLAGS = $(SVNDEF) -O3 -xT 

####all: swift swift_im03 swift_im04 swift_driver swift_testfile swift_window_driver
all: swift 

swift_im03: $(CPPFILES) $(IMAGEANALYSIS0.3) $(HFILES) swift_main.cpp 
	    g++ swift_main.cpp  $(IMAGEANALYSIS0.3) $(CPPFLAGS) -I./ImageAnalysis0.3 $(IFILES) $(LFILES) $(CPPFILES) -o swift_im03

swift_im04: $(CPPFILES) $(IMAGEANALYSIS0.4) $(HFILES) swift_main.cpp 
	    g++ swift_main.cpp  $(IMAGEANALYSIS0.4) $(CPPFLAGS) -I./ImageAnalysis0.4 $(IFILES) $(LFILES) $(CPPFILES) -o swift_im04

swift: $(CPPFILES) $(HFILES) $(SWIFT_IMAGEANALYSIS) swift_main.cpp 
	g++ swift_main.cpp $(CPPFLAGS) -I./SwiftImageAnalysis $(IFILES) $(SWIFT_LFILES) $(CPPFILES) -o $@

swift_intel: $(CPPFILES) $(HFILES) $(SWIFT_IMAGEANALYSIS) swift_main.cpp 
	icpc swift_main.cpp $(INTEL_CPPFLAGS) -I./SwiftImageAnalysis $(IFILES) $(SWIFT_LFILES) $(CPPFILES) -o swift

swift_driver: $(CPPFILES) $(HFILES) $(SWIFT_IMAGEANALYSIS) swift_driver.cpp 
	g++ swift_driver.cpp $(CPPFLAGS) -I./SwiftImageAnalysis $(IFILES) $(SWIFT_LFILES) -o $@

swift_testfile: $(CPPFILES) $(HFILES) $(SWIFT_IMAGEANALYSIS) swift_testfile.cpp 
	g++ swift_testfile.cpp $(CPPFLAGS) -I./SwiftImageAnalysis $(IFILES) $(SWIFT_LFILES) -o $@

swift_window_driver: $(CPPFILES) $(HFILES) $(SWIFT_IMAGEANALYSIS) swift_window_driver.cpp 
	g++ swift_window_driver.cpp $(CPPFLAGS) -I./SwiftImageAnalysis $(IFILES) $(SWIFT_LFILES) -o $@

