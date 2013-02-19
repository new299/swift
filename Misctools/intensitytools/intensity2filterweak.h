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

#ifndef INTENSITYTOOLS_FILTERWEAK_H
#define INTENSITYTOOLS_FILTERWEAK_H

#include "intensityreader.h"
#include "intensitywriter.h"
#include <vector>
#include <iostream>

using namespace std;


void intensity2average(char *input_intensityfilename,vector<double> &a_average,
                                                     vector<double> &g_average,
                                                     vector<double> &c_average,
                                                     vector<double> &t_average) {
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
  
  input_intensity_file.close();
  
  for(vector<vector<double> >::iterator  i= average.begin();i != average.end();i++) {
    for(vector<double>::iterator i1 = (*i).begin();i1 != (*i).end();i1++) {
      (*i1) = (*i1)/count;
    }
  }

  a_average = average[ReadIntensity::base_a];
  t_average = average[ReadIntensity::base_t];
  g_average = average[ReadIntensity::base_g];
  c_average = average[ReadIntensity::base_c];
}

// Extract a single cycle intensity from an intensity file and dump to a new intensity file
// NOTE: Cycles start at 0.
void intensity2filterweak(char *input_intensityfilename,char *input_noisefilename,char *output_intensityfilename,char *output_noisefilename,int cycles_threshold) {

  vector<double> average_a;
  vector<double> average_t;
  vector<double> average_g;
  vector<double> average_c;
  intensity2average(input_intensityfilename,average_a,average_t,average_g,average_c);

  // Intensity file reader
  IntensityReader input_intensity_file (input_intensityfilename);
  IntensityReader input_noise_file     (input_noisefilename);
  IntensityWriter output_intensity_file(output_intensityfilename);
  IntensityWriter output_noise_file    (output_noisefilename);

  for(;!input_intensity_file.eof();) {

    // An intensitysequence is a single line of the intensity file
    IntensitySequence is = input_intensity_file.get_next();
    IntensitySequence in = input_noise_file.get_next();

    if(!input_intensity_file.eof()) {
      bool badread=false;

      int badcycles=0;
      for(int i= 0;i < is.intensities.size();i++) {
        if((is.intensities[i].getbase(ReadIntensity::base_a) < (average_a[i]/2)) &&
           (is.intensities[i].getbase(ReadIntensity::base_t) < (average_t[i]/2)) &&
           (is.intensities[i].getbase(ReadIntensity::base_g) < (average_g[i]/2)) &&
           (is.intensities[i].getbase(ReadIntensity::base_c) < (average_c[i]/2))) {
          badcycles++;
        }
      }
 
      if(badcycles < cycles_threshold) {
        output_intensity_file.write(is);
        output_noise_file.write(in);
      }
    }
  }

  input_intensity_file.close();
  output_intensity_file.close();

  input_noise_file.close();
  output_noise_file.close();
}


#endif
