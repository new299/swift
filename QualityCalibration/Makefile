
SVNDEF := -D'SVN_REV="$(shell svnversion -n .)"'

all: generate_table quality_remap

generate_table: generate_table.cpp *.h
	g++ -O3 $(SVNDEF) generate_table.cpp -I. -o generate_table

quality_remap: quality_remap.cpp *.h
	g++ -O3 $(SVNDEF) quality_remap.cpp -I. -o quality_remap
