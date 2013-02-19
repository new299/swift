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

#ifndef INTENSITYTOOLS_INTENSITY2CORRECTTOCYCLE2_H
#define INTENSITYTOOLS_INTENSITY2CORRECTTOCYCLE2_H

#include "intensityreader.h"
#include "intensitywriter.h"
#include <vector>
#include <iostream>

using namespace std;


// Extract a single cycle intensity from an intensity file and dump to a new intensity file
// NOTE: Cycles start at 0.
void intensity2correcttocycle2(char *input_intensityfilename,char *output_intensityfilename) {

  // Intensity file reader
  IntensityReader input_intensity_file(input_intensityfilename);
  IntensityWriter output_intensity_file(output_intensityfilename);

  int count=0;
  
  
  IntensitySequence is = input_intensity_file.get_next();
  int readsize = is.intensities.size();

  vector<vector <double> > average(readsize,vector<double>(4,0));
  vector<vector <double> > max(readsize,vector<double>(4,0));
  vector<vector <double> > min(readsize,vector<double>(4,0));
  
  vector<vector <bool> > first(readsize,vector<bool>(4,true));
  
  bool first_outer=true;
  for(;!input_intensity_file.eof();) {
   
    // An intensitysequence is a single line of the intensity file
    if(!first_outer) is = input_intensity_file.get_next();
      else first_outer=false;

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
 
  cout << "CORRECTION" << endl;
  for(int b = 0;b < ReadIntensity::base_count;b++) {
    for(int n=0;n<readsize;n++) {
      cout << average[n][b]-average[1][b] << " ";
    }
    cout << endl;
  }

  input_intensity_file.reopen();
  for(;!input_intensity_file.eof();) {
   
    // An intensitysequence is a single line of the intensity file
    IntensitySequence is = input_intensity_file.get_next();
   
    if(!input_intensity_file.eof()) {
      for(int i = 0;i < is.intensities.size();i++) {
        for(ReadIntensity::Base_Type n=0;n<ReadIntensity::base_count;n++) {
          double base = is.intensities[i].getbase(n);
          double origbase = base;
          base = base-(average[i][n]-average[1][n]);
          is.intensities[i].setbase(n,base);
        }
      }
      output_intensity_file.write(is);
    }
  }

  input_intensity_file.close();
  output_intensity_file.close();
}

#endif
