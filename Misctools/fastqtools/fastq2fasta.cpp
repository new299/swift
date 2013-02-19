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

using namespace std;

int main(int argc,char **argv) {

  if(argc < 3) {
    cout << "fastqreader <fastq file> <fasta file>" << endl;
   return 0;
  }

  FastqReader fastq_in(argv[1]); 

  fastq_in.open();
  ofstream outfile(argv[2]);

  cout << "Input file  : " << argv[1] << endl;
  cout << "Output file : " << argv[2] << endl;

  bool eof=false;
  
  for(;!eof;) {
    ScoredSequence s = fastq_in.next_line(eof);
   
    if(!eof) {
      outfile << ">" << s.id << endl;

      outfile << s.sequence;

      outfile << endl;
    }
  }

  fastq_in.close();
  outfile.close();
}
