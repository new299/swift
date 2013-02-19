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

#include "fastareader.h"
#include "fastqwriter.h"
#include "scoredsequence.h"

#include <string>
#include <iostream>
#include <fstream>
#include "stringify.h"

using namespace std;

int main(int argc,char **argv) {

  if(argc < 4) {
    cout << "fasta2fastq <input> <seq output> <quality score>" << endl;
    return 0;
  }

  FastaReader fasta_in(argv[1]); 
  FastqWriter fastq_out(argv[2]);
  
  fasta_in.open();
  fastq_out.open();

  bool eof=false;
  
  for(;!eof;) {
    Sequence s = fasta_in.next_sequence(eof);
    ScoredSequence ss;
    ss.sequence = s.sequence;
    ss.id       = s.id;
    ss.set_quality(convertTo<int>(argv[3]));
    
    fastq_out.write(ss);
  }

  fasta_in.close();
  fastq_out.close();
}
