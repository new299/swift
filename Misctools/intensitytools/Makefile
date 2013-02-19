RM = /bin/rm -f

C++ = g++
CCFLAGS = -g -O3

LDFLAGS = 

all: intensity2pairwiseplot intensity2cycle intensity2purityfilter intensity2normalise intensity2normaliseavg intensity2correcttocycle2 intensity2average intensity2removezero intensity2xybins intensity2flowcellremove intensity2filterweak intensity2filterweakstrong intensity2clean
	

intensity2pairwiseplot: intensity2pairwiseplot.h intensity2pairwiseplot.cpp intensityreader.h intensity.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2pairwiseplot.cpp -o intensity2pairwiseplot $(LDFLAGS)

ipar2gapipelineint: ipar2gapipelineint.cpp iparintensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) ipar2gapipelineint.cpp -o ipar2gapipelineint $(LDFLAGS)

intensity2cycle: intensity2cycle.cpp intensity2cycle.h intensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2cycle.cpp -o intensity2cycle $(LDFLAGS)

intensity2addconst: intensity2addconst.cpp intensityfilter_constant.h intensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2addconst.cpp -o intensity2addconst $(LDFLAGS)

intensity2purityfilter: intensity2purityfilter.cpp intensityfilter_purity.h intensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2purityfilter.cpp -o intensity2purityfilter $(LDFLAGS)

intensity2xybins: intensity2xybins.cpp intensity2xybins.h intensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2xybins.cpp -o intensity2xybins $(LDFLAGS)

intensity2clean: intensity2clean.cpp intensityfilter_flowcell.h intensityfilter_negativezero.h intensity2xybins.h intensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2clean.cpp -o intensity2clean $(LDFLAGS)

intensity2flowcellremove: intensity2flowcellremove.cpp intensityfilter_flowcell.h intensity2xybins.h intensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2flowcellremove.cpp -o intensity2flowcellremove $(LDFLAGS)

intensity2filterblobs: intensity2filterblobs.cpp intensity2filterblobs.h intensityreader.h intensitywriter.h intensity.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2filterblobs.cpp -o intensity2filterblobs $(LDFLAGS)

intensity2removezero: intensity2removezero.cpp intensity2removezero.h intensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2removezero.cpp -o intensity2removezero $(LDFLAGS)

intensity2normalise: intensity2normalise.cpp intensity2normalise.h intensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2normalise.cpp -o intensity2normalise $(LDFLAGS)

intensity2normaliseavg: intensity2normaliseavg.cpp intensity2normaliseavg.h intensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2normaliseavg.cpp -o intensity2normaliseavg $(LDFLAGS)

intensity2correcttocycle2: intensity2correcttocycle2.cpp intensity2correcttocycle2.h intensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2correcttocycle2.cpp -o intensity2correcttocycle2 $(LDFLAGS)

intensity2average: intensity2average.cpp intensity2average.h intensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2average.cpp -o intensity2average $(LDFLAGS)

intensity2filterweak: intensity2filterweak.cpp intensity2filterweak.h intensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2filterweak.cpp -o intensity2filterweak $(LDFLAGS)

intensity2filterweakstrong: intensity2filterweakstrong.cpp intensity2filterweakstrong.h intensityreader.h intensitywriter.h intensity.h 
	@echo --- $@ ---
	$(C++) $(CCFLAGS) intensity2filterweakstrong.cpp -o intensity2filterweakstrong $(LDFLAGS)


install: intensity2pairwiseplot intensity2cycle intensity2purityfilter
	@echo --- $@ ---
	cp intensity2pairwiseplot ../bin
	cp intensity2cycle ../bin
	cp intensity2purityfilter ../bin

clean:
	@echo --- $@ ---
	-$(RM) *.o core
