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

class FastqWriter {
public:

  string filename;
  ofstream fastq_handle;
  int quality_conversion; // This value is subtracted for the ascii character value to obtain the quality score.
  bool textmode;

  FastqWriter() {
  }

  FastqWriter(string filename_in) : filename(filename_in), textmode(false) {
    quality_conversion = 33;
  }

  void set_textmode() {
    textmode = true;
  }

  void qualityoffset(int qualityoffset_in) {
    quality_conversion = qualityoffset_in;
  }


  void write(const ScoredSequence &s) {
    fastq_handle << "@" << s.id << endl;
    fastq_handle << s.sequence << endl;

    fastq_handle << "+" << s.id << endl;
    if(!textmode) fastq_handle << quality_int_to_ascii(s.quality) << endl;
    else {
      for(ScoredSequence::quality_type::const_iterator i = s.quality.begin();i != s.quality.end();i++) {
        fastq_handle << *i << " ";
      }
      fastq_handle << endl;
    }
  }

  void open(string filename_in) {
    filename = filename_in;
    textmode = false;
    open();
  }

  void open() {
    fastq_handle.open(filename.c_str());
  }

  void close() {
    fastq_handle.close();
  }

  string quality_int_to_ascii(const vector<int> qual_in) {
    string s;
    for(vector<int>::const_iterator i = qual_in.begin();i != qual_in.end();i++) {
      s.push_back((char) (*i)+quality_conversion);
    }

    return s;
  }

};
