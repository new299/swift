/*
    Swift (c) 2008 Genome Research Ltd.
    Authors: Nava Whiteford and Tom Skelly (new@sgenomics.org ts6@sanger.ac.uk)

    This file is part of Swift (http://swiftng.sourceforge.net).

    Swift is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Swift is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Swift.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "fastqreader.h"
#include "fastqwriter.h"

#include <string>

#include <stdio.h>
#include <stdlib.h>

using namespace std;

int main(int argc,char **argv) {

  if(argc < 3) {
    cout << "fastqascii2txtqual <fastq file> <fastq txt qual file>" << endl;
    cout << "Unpacks quality values" << endl;
    cout << "Optional [quality offset] argument, allows you to select the offset from the ascii value (i.e. ASCII-33 = quality value)" << endl;
    cout << "the default is 33 (Phred uses 33, Solexa use 64)." << endl;
    return 0;
  }

  FastqReader fastq_in(argv[1]); 
  FastqWriter fastq_out(argv[2]);

  int qualityoffset = -1;
  if(argc > 3) qualityoffset = atoi(argv[3]);

  if(qualityoffset != -1) {
    fastq_in.qualityoffset(qualityoffset);
    fastq_out.qualityoffset(qualityoffset);
  }

  fastq_in.open();
  fastq_out.open();

  cout << "Input file  : " << argv[1] << endl;
  cout << "Output file : " << argv[2] << endl;

  bool eof=false;
 
  fastq_out.set_textmode();

  for(;!eof;) {
    ScoredSequence s = fastq_in.next_line(eof);
    
    fastq_out.write(s);
  }

  fastq_in.close();
  fastq_out.close();
}
