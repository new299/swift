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

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "scoredsequence.h"

using namespace std;

class FastqReader {
public:

  string filename;
  ifstream fastq_handle;
  int quality_conversion; // This value is subtracted for the ascii character value to obtain the quality score.

  FastqReader(string filename_in) : filename(filename_in) {
    quality_conversion = 33;
  }

  void qualityoffset(int qualityoffset_in) {
    quality_conversion = qualityoffset_in;
  }

  vector<ScoredSequence> parse_all() {
    vector<ScoredSequence> v;
    
    open();

    bool eof = false;
    for(;!eof;) {
      ScoredSequence s = next_line(eof);
      if(eof == false) v.push_back(s);
    }

    close();

    return v;
  }

  void open() {
    fastq_handle.open(filename.c_str());
  }

  // reads a line in to a ScoredSequence
  ScoredSequence next_line(bool &eof) {
    string line;
    ScoredSequence s;
    
    // ID line (strip leading character)
    getline(fastq_handle,line);
    if(fastq_handle.eof()) {eof = true; return s;}

    s.id = line.substr(1);
    
    // Sequence line
    getline(fastq_handle,line);
    if(fastq_handle.eof()) {eof = true; return s;}
    s.sequence = line;
    
    // Another ID line, which I will ignore
    getline(fastq_handle,line);
    if(fastq_handle.eof()) {eof = true; return s;}

    // Quality line
    getline(fastq_handle,line);
    if(fastq_handle.eof()) {eof = true; return s;}
   
    for(int n=0;n<line.length();n++) {
      s.quality.push_back((int) line[n]-quality_conversion);
    }
  }

  void close() {
    fastq_handle.close();
  }

};
