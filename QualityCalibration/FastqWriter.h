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

#ifndef SWIFT_FASTQWRITER_H
#define SEIFT_FASTQWRITER_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "ProbabilitySequence.h"

using namespace std;

template<class base_type>
class FastqWriter {
public:

  string filename;
  ofstream fastq_handle;
  int quality_conversion; // This value is subtracted for the ascii character value to obtain the quality score.
  bool textmode;

  FastqWriter(string filename_in) : filename(filename_in), textmode(false) {
    quality_conversion = 33;
    openok=false;
    open();
  }

  void set_textmode() {
    textmode = true;
  }

  void qualityoffset(int qualityoffset_in) {
    quality_conversion = qualityoffset_in;
  }
 
  int get_quality_offset() {
    return quality_conversion;
  }

  bool is_open() {
    return openok;
  }

  template<class _prec>
  string quality_int_to_ascii(const vector<_prec> qual_in) {
    string s;
    for(typename vector<_prec>::const_iterator i = qual_in.begin();i != qual_in.end();i++) {
      s.push_back(static_cast<char>((*i)+quality_conversion+0.5));
    }

    return s;
  }

  void write(const ProbabilitySequence<base_type> &s) {
    fastq_handle << "@" << s.get_id() << endl;
    fastq_handle << s.get_sequence_string() << endl;

    fastq_handle << "+" << s.get_id() << endl;
    if(!textmode) fastq_handle << quality_int_to_ascii(s.phred_quality()) << endl;
    else {
      for(typename vector<int>::const_iterator i = s.phred_quality().begin();i != s.phred_quality().end();i++) {
        fastq_handle << *i << " ";
      }
      fastq_handle << endl;
    }
  }

  void write(const string &id,
             const string &sequence,
             const string &quality) {
    fastq_handle << "@" << id << endl;
    fastq_handle << sequence << endl;

    fastq_handle << "+" << id << endl;
    fastq_handle << quality << endl;
  }

  void open() {
    fastq_handle.open(filename.c_str());
    openok = fastq_handle.good();
  }

  void close() {
    fastq_handle.close();
  }

  bool openok;

};

#endif
