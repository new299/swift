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

#ifndef SWIFT_FASTPRBREADER_H
#define SWIFT_FASTPRBREADER_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "ProbabilitySequence.h"
#include "stringify.h"
#include <sstream>

using namespace std;

template<class base_type>
class FastPrbReader {
public:

  string filename;
  ifstream          fastprb_handle;
  unsigned int      encode_mode;
  int               quality_offset;
  size_t            probability_count;

  static const unsigned int      encodemode_ascii    = 0;
  static const unsigned int      encodemode_text     = 1;
  static const unsigned int      encodemode_unknown  = 2;

  FastPrbReader(string filename_in) : filename(filename_in) {
    open();
    quality_offset=0;
  }

  void open() {
    fastprb_handle.open(filename.c_str());
    parse_header();
    getline(fastprb_handle,line);
  }

  void close() {
    fastprb_handle.close();
  }

  bool eof() {
    return fastprb_handle.eof();
  }

  void parse_header() {
    fastprb_handle.seekg(0,ios::beg);
    
    ProbabilitySpecification &spec = base_type::get_probability_specification();
    
    probability_count=0;
    bool stop=false;
    for(;!stop;) {
      string line;
      int beforeline = fastprb_handle.tellg();
      getline(fastprb_handle,line);

      if(line.compare(0,8,string("\\Version")) == 0) {
        
        string::size_type pos = line.find_first_of("=");
        string version_spec = line.substr(pos+1);
        
        // parse version string
        pos = version_spec.find_first_of(":");
        string id_str = version_spec.substr(0,pos+1);

        if(id_str.compare(0,5,string("FASTPRB")) == 0) {
          cerr << "ID String OK" << endl;
        } else cerr << "ID String BAD" << endl;
        
        version_spec = version_spec.substr(pos+1);

        // parse version string
        pos = version_spec.find_first_of(":");

        string version_number_str = version_spec.substr(0,pos);
        cout << "version string: " << version_number_str << endl;
        int version = convertTo<int>(version_number_str);
        
        version_spec = version_spec.substr(pos+1);
        
        // set text mode
        pos = version_spec.find_first_of(":");

        if(pos == string::npos) {
          pos = version_spec.size();
        }

        string texttype_str  = version_spec.substr(0,pos+1);

        if(texttype_str.compare(0,1,string("T")) == 0) {
          
          encode_mode = encodemode_text;
          cerr << "textmode" << endl;
        } else
        if(texttype_str.compare(0,1,string("A")) == 0) {
          encode_mode = encodemode_ascii;
          cerr << "ascii encoded mode" << endl;
        } else encode_mode = encodemode_unknown;
      }

      if(line.compare(0,5,string("\\Base")) == 0) {
        
        string::size_type pos = line.find_first_of("=");
        string base_spec = line.substr(pos+1);

        // parse Base spec

        // to next "," grab as int, index number
        pos = base_spec.find_first_of(",");
        string index_str = base_spec.substr(0,pos);

        // cout << "index string: " << index_str << endl;
        int index = convertTo<int>(index_str);
        
        base_spec = base_spec.substr(pos+1);

        // to next "," grab as string, short name
        pos = base_spec.find_first_of(",");
        string shortname = base_spec.substr(0,pos);
        // cout << "shortname string: " << shortname << endl;
        base_spec = base_spec.substr(pos+1);

        // to next "," grab as string, long description to end of line.
        string longname = base_spec;
        // cout << "longname string: " << longname << endl;
        spec.add_specification(index,shortname,longname);
        probability_count++;
      }

      if(line.compare(0,1,string("@")) == 0) {
        fastprb_handle.seekg(beforeline);
        stop = true;
      }
    }
  }

  size_t get_probabilies_count() {
    // Move to start of file
    fastprb_handle.seekg(0,ios::beg);

    // Count lines that start with \Base until we reach at @
    size_t count=0;
    bool stop=false;
    for(;!stop;) {
      string line;
      getline(fastprb_handle,line);

      if(line.compare(0,5,string("\\Base")) == 0) {count++;  }
      if(line.compare(0,1,string("@")     ) == 0) {stop=true;}
    }

    return count;
  }
  
  typename base_type::probability_type get_probability(string &prbstr) {
  
    if(encode_mode == encodemode_ascii) {
      typename base_type::probability_type prbchar = prbstr[0];
      prbstr = prbstr.substr(1);

      prbchar -= quality_offset;

      // magic to convert ascii value in to probability
      prbchar = prbchar/10;
      prbchar = pow(10,prbchar);
      prbchar = prbchar/(prbchar+1);

      return prbchar;

    } else 
    if(encode_mode == encodemode_text) {

      istringstream ss(prbstr);

      string prb;
      ss >> prb;

      if(prbstr.size() > (prb.size()+1)) { prbstr = prbstr.substr(prb.size()+1);}
                                    else { prbstr = ""; }

      return convertTo<typename base_type::probability_type>(prb);
    }
  }

  ProbabilitySequence<base_type> get_sequence() {
    
    ProbabilitySequence<base_type> s;

    string tag;
    
    if(line.compare(0,1,string("@")) == 0) {
      // Grab base tag
      tag = line.substr(1);

      s.set_id(tag);

      // Grab sequence probabilities
      vector<string> probability_strings;
      for(size_t n=0;n< probability_count;n++) {
        getline(fastprb_handle,line);
        probability_strings.push_back(line);
      }

      // Grab first set of probabilities from strings
      for(;probability_strings[0].size() != 0;) {
        base_type prb(probability_count);
        for(size_t n=0;n<probability_count;n++) {
          prb.set_base(n,get_probability(probability_strings[n]));
        }

        s.sequence().push_back(base_type(prb));
      }
    }
    
    getline(fastprb_handle,line);

    return s;
  }
  
  string line;

};

#endif
