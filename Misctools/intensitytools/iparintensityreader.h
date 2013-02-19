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

#ifndef SIMPLECALLER_IPARINTENSITYREADER
#define SIMPLECALLER_IPARINTENSITYREADER

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "intensity.h"

using namespace std;

class IparIntensityReader {

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

  IparIntensityReader() {
  }

  IparIntensityReader(string intensity_filename_in,
                      string positions_filename_in,
		      int lane_in,
		      int tile_in) 
		    : intensity_filename(intensity_filename_in),
		      positions_filename(positions_filename_in),
                      lane(lane_in),
                      tile(tile_in),
                      intensities_idx(0) {

    // We have to load all data in to memory currently, due to the way the information is held in the intensity files.
    load_all_data();

  }

  void load_all_data() {
    ifstream intensity_file(intensity_filename.c_str());
    ifstream positions_file(positions_filename.c_str());

    // 1. Create new intensity sequences at correct positions, index by position in positions file.
    for(size_t n=0;!positions_file.eof();n++) {

      string current_line;
      getline(positions_file,current_line);

      istringstream current_line_ss(current_line);

      double x_position;
      double y_position;

      current_line_ss >> x_position;
      current_line_ss >> y_position;

      IntensitySequence s(static_cast<int>(x_position),static_cast<int>(y_position),lane,tile);

      intensities.push_back(s);
    }

    // 2. Add intensity data to sequences
    size_t cluster_idx = 0;
    for(size_t n=0;!intensity_file.eof();n++) {
      string current_line;
      getline(intensity_file,current_line);

      if(current_line[0] == '#') {
        // end of cycle...
        cout << "hit tag: " << current_line << endl;
        cout << "final cluster idx: " << cluster_idx << endl;
	cluster_idx = 0;
      } else {

        double a_intensity;
        double c_intensity;
        double g_intensity;
        double t_intensity;

	istringstream current_line_ss(current_line);

        current_line_ss >> a_intensity;
	current_line_ss >> c_intensity;
	current_line_ss >> g_intensity;
	current_line_ss >> t_intensity;

        ReadIntensity i(a_intensity,c_intensity,g_intensity,t_intensity);

        intensities[cluster_idx].intensities.push_back(i);

        cluster_idx++;
      }
    }
  }

  IntensitySequence get_next() {
    intensities_idx++;
    return intensities[intensities_idx-1];
   
  }

  bool eof() {
    if(intensities_idx >= intensities.size()) return true;
                                         else return false;
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
 
  size_t intensities_idx;
  vector<IntensitySequence> intensities; 
  string input_filename;
  ifstream input_file;
  string intensity_filename;
  string positions_filename;
  int lane;
  int tile;
};

#endif
