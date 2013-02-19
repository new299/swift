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


#ifndef FASTATOOLS_FASTAREADER_H
#define FASTATOOLS_FASTAREADER_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "sequence.h"

using namespace std;

class FastaReader {
public:

  string filename;
  ifstream fasta_handle;

  string next_contig_name;

  FastaReader(string filename_in) : filename(filename_in) {
  }

  void open() {
    fasta_handle.open(filename.c_str());
    string line;

    // Read inital contig name
    getline(fasta_handle,line);
    if(!fasta_handle.eof()) next_contig_name = line.substr(1);
  }

  // reads a contig in to a Sequence
  Sequence next_sequence(bool &eof) {
    string line;
    Sequence s;
    
    s.id = next_contig_name;
    next_contig_name = "";

    for(;;) {
      // Get a line
      getline(fasta_handle,line);
      if(fasta_handle.eof()) {eof = true; return s;}
    
      // if the line starts with a > then we have reached the end of the contig.
      // set the next contig name to the remainder of this line and return.
      if(line[0] == '>') {
        next_contig_name = line.substr(1);
        return s;
      } else {
        s.sequence += line;
      }
    }
  }

  void close() {
    fasta_handle.close();
  }

};

#endif
