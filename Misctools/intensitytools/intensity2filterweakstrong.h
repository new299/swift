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

#ifndef INTENSITYTOOLS_FILTERWEAKSTRONG_H
#define INTENSITYTOOLS_FILTERWEAKSTRONG_H

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

  vector<vector<int> >       count(4,vector<int>   (is.intensities.size(),0));
  vector<vector <double> > average(4,vector<double>(is.intensities.size(),0));
  
  bool firstloop=true;
  vector<vector <bool> > first(is.intensities.size(),vector<bool>(4,true));
  for(;!input_intensity_file.eof();) {
    if(!firstloop) {
      is = input_intensity_file.get_next();
    } else firstloop=false; 

    // An intensitysequence is a single line of the intensity file
    if(!input_intensity_file.eof()) {
      for(int i = 0;i < is.intensities.size();i++) {
          average[is.intensities[i].max_base()][i] += is.intensities[i].getbase(is.intensities[i].max_base());
          count[is.intensities[i].max_base()][i]++;
      }
    }
  }
  
  input_intensity_file.close();
  
  for(int  i= 0;i < average.size();i++) {
    for(int i1 = 0;i1 < average[i].size();i1++) {
      average[i][i1] = average[i][i1]/count[i][i1];
    }
  }

  a_average = average[ReadIntensity::base_a];
  t_average = average[ReadIntensity::base_t];
  g_average = average[ReadIntensity::base_g];
  c_average = average[ReadIntensity::base_c];
}

// Extract a single cycle intensity from an intensity file and dump to a new intensity file
// NOTE: Cycles start at 0.
void intensity2filterweakstrong(char *input_intensityfilename,char *input_noisefilename,char *output_intensityfilename,char *output_noisefilename,int cycles_threshold) {

  vector<double> average_a;
  vector<double> average_t;
  vector<double> average_g;
  vector<double> average_c;
  intensity2average(input_intensityfilename,average_a,average_t,average_g,average_c);

  cout << "average a: ";
  for(int n=0;n<average_a.size();n++) {
    cout << average_a[n] << " ";
  }
  cout << endl;
  
  cout << "average t: ";
  for(int n=0;n<average_t.size();n++) {
    cout << average_t[n] << " ";
  }
  cout << endl;
  
  cout << "average g: ";
  for(int n=0;n<average_g.size();n++) {
    cout << average_g[n] << " ";
  }
  cout << endl;
  
  cout << "average c: ";
  for(int n=0;n<average_c.size();n++) {
    cout << average_c[n] << " ";
  }
  cout << endl;

  cout << "Cycles Threshold: " << cycles_threshold << endl;

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
        if((is.intensities[i].getbase(ReadIntensity::base_a) < (average_a[i]/3)) &&
           (is.intensities[i].getbase(ReadIntensity::base_t) < (average_t[i]/3)) &&
           (is.intensities[i].getbase(ReadIntensity::base_g) < (average_g[i]/3)) &&
           (is.intensities[i].getbase(ReadIntensity::base_c) < (average_c[i]/3))) {
          badcycles++;
        }
        
        int badint=0;
        if(is.intensities[i].getbase(ReadIntensity::base_a) > (average_a[i]*3)) badint++; 
        if(is.intensities[i].getbase(ReadIntensity::base_t) > (average_t[i]*3)) badint++;
        if(is.intensities[i].getbase(ReadIntensity::base_g) > (average_g[i]*3)) badint++;
        if(is.intensities[i].getbase(ReadIntensity::base_c) > (average_c[i]*3)) badint++;
        if(badint >= 2) badcycles++;
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
