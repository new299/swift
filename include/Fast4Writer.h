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

#ifndef SWIFT_FAST4WRITER_H
#define SEIFT_FAST4WRITER_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "ProbabilitySequence.h"

using namespace std;

template<class _basetype=int, class _prbprec=double>
class Fast4Writer {
public:

  string filename;
  ofstream fast4_handle;
  int quality_conversion; // This value is subtracted for the ascii character value to obtain the quality score.
  bool textmode;

  Fast4Writer(string filename_in) : filename(filename_in), textmode(true) {
    quality_conversion = 74;
  }

  void set_textmode() {
    textmode = true;
  }

  void qualityoffset(int qualityoffset_in) {
    quality_conversion = qualityoffset_in;
  }
 
  string quality_int_to_ascii(const vector<_prbprec> qual_in) {
    string s;
    for(typename vector<_prbprec>::const_iterator i = qual_in.begin();i != qual_in.end();i++) {
      s.push_back(static_cast<char>((*i)+quality_conversion+0.5));
    }

    return s;
  }


  string get_fast4_header() {
    ostringstream s;
    
    s << "\\Version=FAST4:1";
    if(textmode) { s << ":T"; }
            else { s << ":A"; }

    s << endl;

    // Base format
    for(size_t n=0;n<ProbabilitySequence<_basetype,_prbprec>::base_count;n++) {
      s << "\\Base=" << n << "," << ProbabilitySequence<_basetype,_prbprec>::base_name[n] << "," << ProbabilitySequence<_basetype,_prbprec>::base_descriptions[n] << endl;
    }
    
    if(!textmode) s << "\\QualityOffset=" << quality_conversion << endl;

    return s.str();
  }

  string get_fast4_textmode_string_sequence(const ProbabilitySequence<_basetype,_prbprec> &p) const {

    ostringstream s;

    s << "@" << p.get_id() << endl;

    for(size_t n=0;n<ProbabilitySequence<_basetype,_prbprec>::base_count;n++) {
      for(typename ProbabilitySequence<_basetype,_prbprec>::probability_sequence_type::const_iterator i=p.const_sequence().begin();i != p.const_sequence().end();i++) {
        s << fixed;
        double c = (*i)[n];
        if(c == -0) c=0;

        s << setprecision(4) << setw(4) << c << " ";
      }
      s << endl;
    }

    return s.str();
  }
  
  string get_fast4_ascii_string_sequence(const ProbabilitySequence<_basetype,_prbprec> &p) const {

    ostringstream s;

    s << "@" << p.get_id() << endl;

    for(size_t n=0;n<ProbabilitySequence<_basetype,_prbprec>::base_count;n++) {
      for(typename ProbabilitySequence<_basetype,_prbprec>::probability_sequence_type::const_iterator i=p.const_sequence().begin();i != p.const_sequence().end();i++) {
        int c = static_cast<int>(10*(log10(((*i)[n])/(1-(*i)[n]))));
        if(c>= 50) c= 50;
        if(c<=-40) c=-40;

        c += quality_conversion;

        s << static_cast<char>(c);
      }
      s << endl;
    }

    return s.str();
  }

  void write(const ProbabilitySequence<_basetype,_prbprec> &s) {
    if( textmode) fast4_handle << get_fast4_textmode_string_sequence(s);
    if(!textmode) fast4_handle << get_fast4_ascii_string_sequence(s);
  }

  void open() {
    fast4_handle.open(filename.c_str());

    fast4_handle << get_fast4_header();
  }

  void close() {
    fast4_handle.close();
  }
};

#endif
