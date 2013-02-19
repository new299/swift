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

#ifndef SIMPLECALLER_INTENSITYREADER
#define SIMPLECALLER_INTENSITYREADER

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "intensity.h"

using namespace std;

class IntensityReader {

public:
  
  // make private, encapsulate 
  vector<double> max_a;
  vector<double> min_a;
  vector<double> max_t;
  vector<double> min_t;
  vector<double> max_g;
  vector<double> min_g;
  vector<double> max_c;
  vector<double> min_c;
 
  int max_readlength;

  int total;

  IntensityReader() {
  }

  IntensityReader(string filename) : input_filename(filename) {
    input_file.open(filename.c_str());
  }

  IntensitySequence get_next() {
    return read_sequence_intensity();
  }

  bool eof() {
    return input_file.eof();
  }

  void close() {
    input_file.close();
  }

  void reopen() {
    input_file.close();
    input_file.open(input_filename.c_str());
  }

  void calc_stats() {
    max_a.clear();
    max_a.insert(max_a.begin(),100,0);
    min_a.insert(min_a.begin(),100,0);
    max_t.insert(max_t.begin(),100,0);
    min_t.insert(min_t.begin(),100,0);
    max_g.insert(max_g.begin(),100,0);
    min_g.insert(min_g.begin(),100,0);
    max_c.insert(max_c.begin(),100,0);
    min_c.insert(min_c.begin(),100,0);

    total=0;

    int count=0;
    max_readlength = 0;
    bool first=true;
    for(;!eof();) {
      IntensitySequence s = get_next();
   
      int position=0;
      for(IntensitySequence::intensities_type::iterator i=s.intensities.begin();i != s.intensities.end();i++,position++) {

        if(position > max_readlength) max_readlength = position;
	if(first == true) {
	  max_a[position] = (*i).getbase(ReadIntensity::base_a);
	  max_t[position] = (*i).getbase(ReadIntensity::base_t);
	  max_g[position] = (*i).getbase(ReadIntensity::base_g);
	  max_c[position] = (*i).getbase(ReadIntensity::base_c);
	  
	  min_a[position] = (*i).getbase(ReadIntensity::base_a);
	  min_t[position] = (*i).getbase(ReadIntensity::base_t);
	  min_g[position] = (*i).getbase(ReadIntensity::base_g);
	  min_c[position] = (*i).getbase(ReadIntensity::base_c);
	  
	  first=false;
	}
	
	if( (*i).getbase(ReadIntensity::base_a) > max_a[position]) { max_a[position] = (*i).getbase(ReadIntensity::base_a); }
	if( (*i).getbase(ReadIntensity::base_t) > max_t[position]) { max_t[position] = (*i).getbase(ReadIntensity::base_t); }
	if( (*i).getbase(ReadIntensity::base_g) > max_g[position]) { max_g[position] = (*i).getbase(ReadIntensity::base_g); }
	if( (*i).getbase(ReadIntensity::base_c) > max_c[position]) { max_c[position] = (*i).getbase(ReadIntensity::base_c); }
	
	if( (*i).getbase(ReadIntensity::base_a) < min_a[position]) { min_a[position] = (*i).getbase(ReadIntensity::base_a); }
	if( (*i).getbase(ReadIntensity::base_t) < min_t[position]) { min_t[position] = (*i).getbase(ReadIntensity::base_t); }
	if( (*i).getbase(ReadIntensity::base_g) < min_g[position]) { min_g[position] = (*i).getbase(ReadIntensity::base_g); }
	if( (*i).getbase(ReadIntensity::base_c) < min_c[position]) { min_c[position] = (*i).getbase(ReadIntensity::base_c); }

	count++;
      }
    }

    total=count;
  }

private:
  IntensitySequence read_sequence_intensity() {
    IntensitySequence s;
    input_file >> s.lane;
    input_file >> s.tile;
    input_file >> s.x;
    input_file >> s.y;

    string str;
    getline(input_file, str);
    istringstream iss(str);
    
    for(;!iss.eof();) {
      ReadIntensity r;
      r << iss;
      s.intensities.push_back(r);
    }

    return s;
  }
  
  string input_filename;
  ifstream input_file;
};

#endif
