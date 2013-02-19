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

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc,char **argv) {

  if(argc < 3) {
    cout << "fasta2flat <input> <seq output> <contig names output> <contig positions output>" << endl;
    return 0;
  }

  FastaReader fasta_in(argv[1]); 
  ofstream    flat_out(argv[2]);
  ofstream    contignames_out(argv[3]);
  ofstream    contigpositions_out(argv[4]);
 
  fasta_in.open();
  

  bool eof=false;
  
  int position=0;
  for(;!eof;) {
    Sequence s = fasta_in.next_sequence(eof);
    position += s.sequence.length();

    flat_out        << s.sequence << endl;
    contignames_out << s.id << endl;
    contigpositions_out << position << endl;
  }

  fasta_in.close();
  flat_out.close();
  contignames_out.close();
  contigpositions_out.close();
}
