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

#ifndef INTENSITYTOOLS_FILTERBLOBS_H
#define INTENSITYTOOLS_FILTERBLOBS_H

#include "intensityreader.h"
#include "intensitywriter.h"
#include <vector>
#include <iostream>

using namespace std;


// Extract a single cycle intensity from an intensity file and dump to a new intensity file
// NOTE: Cycles start at 0.
void intensity2filterblobs(char *input_intensityfilename,char *input_noisefilename,char *output_intensityfilename,char *output_noisefilename,int cycles_threshold) {

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

      double old_max_intensity=0;
      ReadIntensity::Base_Type old_max_base     = ReadIntensity::base_a;
      
      double current_max_intensity = is.intensities[0].max_intensity();
      ReadIntensity::Base_Type current_max_base = is.intensities[0].max_base();
      int badcycles=0;

      for(int i= 1;i < is.intensities.size();i++) {
        old_max_intensity = current_max_intensity;
        old_max_base = current_max_base;
        
        current_max_intensity = is.intensities[i].max_intensity();
        current_max_base      = is.intensities[i].max_base();

        if((current_max_base == old_max_base) && (current_max_intensity >= old_max_intensity)) badcycles++;
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
