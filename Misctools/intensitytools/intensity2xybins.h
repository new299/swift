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

#ifndef INTENSITYTOOLS_INTENSITY2XYBINS_H
#define INTENSITYTOOLS_INTENSITY2XYBINS_H

#include "intensityreader.h"
#include "intensitywriter.h"
#include <vector>
#include <iostream>

using namespace std;

void intensity2xybins(IntensityReader &input_intensity_file,int bin_count,double &bin_size, vector<int> &xbins,vector<int> &ybins) {

  // find max of x and y (we'll make an equal number of bins in each direction)
  bool first=true;
  double max=0;
  for(;!input_intensity_file.eof();) {
   
    // An intensitysequence is a single line of the intensity file
    IntensitySequence is = input_intensity_file.get_next();

    if(!input_intensity_file.eof()) {
      if(is.x > max || first) {max = is.x; first=false;}
      if(is.y > max) max = is.y;
    }
  }

  input_intensity_file.reopen();

  // count strips across tile
  xbins.clear();
  ybins.clear();
  xbins.insert(xbins.begin(),bin_count,0);
  ybins.insert(ybins.begin(),bin_count,0);
  bin_size = (max/bin_count);
  for(;!input_intensity_file.eof();) {
   
    // An intensitysequence is a single line of the intensity file
    IntensitySequence is = input_intensity_file.get_next();

    if(!input_intensity_file.eof()) {
      // a little messy, but basically find the bin position (rounded) and increment
      xbins[static_cast<int>((is.x/(max/bin_count))+0.5)]++;
      ybins[static_cast<int>((is.y/(max/bin_count))+0.5)]++;
    }
  }
  
  input_intensity_file.close();
}

#endif
