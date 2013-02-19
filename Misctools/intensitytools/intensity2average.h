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

#ifndef INTENSITYTOOLS_AVERAGE_H
#define INTENSITYTOOLS_AVERAGE_H

#include "intensityreader.h"
#include "intensitywriter.h"
#include <vector>
#include <iostream>

using namespace std;

// Extract a single cycle intensity from an intensity file and dump to a new intensity file
// NOTE: Cycles start at 0.
void intensity2average(char *input_intensityfilename,char *output_filename) {

  // Intensity file reader
  IntensityReader input_intensity_file(input_intensityfilename);
  IntensitySequence is = input_intensity_file.get_next();

  int count=0;
  vector<vector <double> > average(is.intensities.size(),vector<double>(4,0));
  vector<vector <double> > max(is.intensities.size(),vector<double>(4,0));
  vector<vector <double> > min(is.intensities.size(),vector<double>(4,0));
  
  bool firstloop=true;
  vector<vector <bool> > first(is.intensities.size(),vector<bool>(4,true));
  for(;!input_intensity_file.eof();) {
    if(!firstloop) {
      is = input_intensity_file.get_next();
    } else firstloop=false; 

    // An intensitysequence is a single line of the intensity file
    if(!input_intensity_file.eof()) {
      for(int i = 0;i < is.intensities.size();i++) {
        for(ReadIntensity::Base_Type n=0;n<ReadIntensity::base_count;n++) {
          average[i][n] += is.intensities[i].getbase(n);
          if((is.intensities[i].getbase(n) > max[i][n]) || first[i][n]) {max[i][n] = is.intensities[i].getbase(n); first[i][n]=false;}
          if((is.intensities[i].getbase(n) < min[i][n]) || first[i][n]) {min[i][n] = is.intensities[i].getbase(n); first[i][n]=false;}
        }
      }
      count++;
    }
  }
  
  for(vector<vector<double> >::iterator  i= average.begin();i != average.end();i++) {
    for(vector<double>::iterator i1 = (*i).begin();i1 != (*i).end();i1++) {
      (*i1) = (*i1)/count;
    }
  }

  ofstream a_average_file((string(output_filename) + "_a").c_str());
  ofstream t_average_file((string(output_filename) + "_t").c_str());
  ofstream g_average_file((string(output_filename) + "_g").c_str());
  ofstream c_average_file((string(output_filename) + "_c").c_str());

  for(int i=0;i < average.size();i++) {
    for(int b=0;b < average[i].size();b++) {
      if(b == ReadIntensity::base_a) a_average_file << i << " " << average[i][b] << endl;
      if(b == ReadIntensity::base_t) t_average_file << i << " " << average[i][b] << endl;
      if(b == ReadIntensity::base_g) g_average_file << i << " " << average[i][b] << endl;
      if(b == ReadIntensity::base_c) c_average_file << i << " " << average[i][b] << endl;
    }
  }

  a_average_file.close();
  t_average_file.close();
  g_average_file.close();
  c_average_file.close();
  input_intensity_file.close();
}

#endif
