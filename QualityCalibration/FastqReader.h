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

#ifndef SWIFT_FASTQREADER_H
#define SWIFT_FASTQREADER_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "ProbabilitySequence.h"
#include "stringify.h"
#include <sstream>

using namespace std;

template<class base_type>
class FastqReader {
public:

  string filename;
  string last_quality;
  ifstream          fasta_handle;
  static const int  fasta_probabilities = 5;
  bool openok;

  int qualityoffset;
  string next_contig_name;

  FastqReader(string filename_in) : filename(filename_in) {
    
    ProbabilitySpecification &spec = base_type::get_probability_specification();

    qualityoffset = 33;
    spec.clear();
    spec.add_specification(0,"A","Adenine" );
    spec.add_specification(1,"C","Cytosine");
    spec.add_specification(2,"G","Guanine" );
    spec.add_specification(3,"T","Thymine" );

    open();
  }

  void close() {
    fasta_handle.close();
  
  }

  const string &get_last_quality_string() {
    return last_quality;
  }

  const int get_quality_offset() {
    return qualityoffset;
  }

  bool open() {
    openok=false;
    fasta_handle.open(filename.c_str());
    if(fasta_handle.good()) openok=true;
    string line;

    // Read initial contig name
    if(openok) {
      getline(fasta_handle,line);
      if(!fasta_handle.eof()) next_contig_name = line.substr(1);
    }

    return openok;
  }

  bool is_open() {
    return openok;
  }

  // reads a contig in to a Sequence
  ProbabilitySequence<base_type> get_sequence() {

    ProbabilitySequence<base_type> sprb;

    string line;
    string s;
    string q;
    
    sprb.set_id(next_contig_name);
    next_contig_name = "";

    bool first=true;

    int linecount=0;
    for(;;) {
      // Get a line
      getline(fasta_handle,line);
      first=false;

      // if the line starts with a > then we have reached the end of the contig.
      // set the next contig name to the remainder of this line and return.
      //
      if(linecount % 4 == 1) {}

      if(linecount % 4 == 3 || (fasta_handle.eof()) ) {
        if(!fasta_handle.eof()) next_contig_name = line.substr(1);

        // translate s in to a ProbabilitySequence<>
        for(int n=0;n<s.length();n++) {
          //cout << "s: " << s << endl;
          //cout << "q: " << q << endl;
          string shortname = s.substr(n,1);
          string shortprb  = q.substr(n,1);
          size_t idx = base_type::get_probability_specification().get_index(shortname);
          typename base_type::probability_type p = asciiprb_to_double(shortprb);

          base_type b(fasta_probabilities);
          for(size_t i=0;i<base_type::get_probability_specification().get_probability_count();i++) {
            if(i == idx) { b.set_base(i,p); }
                    else { b.set_base(i,(1-p)/fasta_probabilities); }
          }
          sprb.sequence().push_back(b);
        }

        //cout << "returning" << endl;
        return sprb;
      } 

      if(linecount % 4 == 0) {
        s = line;
      }

      if(linecount % 4 == 2) {
        q = line;
        last_quality = q;
      }

      linecount++;
    }
  }

  double asciiprb_to_double(string q) {
    double prbchar = q[0];

    prbchar = prbchar - qualityoffset;
    
    // magic to convert ascii value in to probability
    prbchar = prbchar/10;
    prbchar = pow(10,prbchar);
    prbchar = prbchar/(prbchar+1);

    return prbchar;  
  }

  bool eof() {
    return fasta_handle.eof();
  }
};

#endif
