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

  if(argc < 2) {
    cout << "fastacount <input>" << endl;
    return 0;
  }

  FastaReader fasta_in(argv[1]); 
 
  fasta_in.open();
  

  bool eof=false;
  
  int a_count=0;
  int t_count=0;
  int g_count=0;
  int c_count=0;
  int n_count=0;

  int position=0;
  for(;!eof;) {
    Sequence s = fasta_in.next_sequence(eof);

    for(int n=0;n<s.sequence.length();n++) {
      if(s.sequence[n] == 'A') a_count++; else
      if(s.sequence[n] == 'T') t_count++; else
      if(s.sequence[n] == 'G') g_count++; else
      if(s.sequence[n] == 'C') c_count++; else
        n_count++;
    }
  }

  cout << "RAW" << endl;
  cout << "A Count: " << a_count << endl;
  cout << "T Count: " << t_count << endl;
  cout << "G Count: " << g_count << endl;
  cout << "C Count: " << c_count << endl;
  cout << "N Count: " << n_count << endl;
  cout << endl;

  int total = a_count+t_count+g_count+c_count+n_count;

  cout << "Normalised" << endl;
  cout << "A Count: " << static_cast<float>(a_count)/static_cast<float>(total) << endl;
  cout << "T Count: " << static_cast<float>(t_count)/static_cast<float>(total) << endl;
  cout << "G Count: " << static_cast<float>(g_count)/static_cast<float>(total) << endl;
  cout << "C Count: " << static_cast<float>(c_count)/static_cast<float>(total) << endl;
  cout << "N Count: " << static_cast<float>(n_count)/static_cast<float>(total) << endl;

  fasta_in.close();
}
