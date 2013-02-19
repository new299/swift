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

#ifndef INTENSITYTOOLS_INTENSITY2CYCLE_H
#define INTENSITYTOOLS_INTENSITY2CYCLE_H

#include "intensityreader.h"
#include "intensitywriter.h"
#include <vector>
#include <iostream>

using namespace std;


// Extract a single cycle intensity from an intensity file and dump to a new intensity file
// NOTE: Cycles start at 0.
void intensity2cycle(char *input_intensityfile,char *output_intensityfile,int cycle) {

  // Intensity file reader
  IntensityReader input_intensity_file(input_intensityfile);
  IntensityWriter output_intensity_file(output_intensityfile);

  for(;!input_intensity_file.eof();) {
   
    // An intensitysequence is a single line of the intensity file
    IntensitySequence is = input_intensity_file.get_next();
    
    if(!input_intensity_file.eof()) {
      if(is.intensities.size() > cycle) {
        IntensitySequence is_new(is.x,is.y,is.lane,is.tile);
	is_new.intensities.push_back(is.intensities[cycle]);
        output_intensity_file.write(is_new);
      }
    }
  }

  input_intensity_file.close();
  output_intensity_file.close();
}

#endif
