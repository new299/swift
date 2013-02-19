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

#ifndef INTENSITYTOOLS_INTENSITY2REMOVEZERO_H
#define INTENSITYTOOLS_INTENSITY2REMOVEZERO_H

#include "intensityreader.h"
#include "intensitywriter.h"
#include <vector>
#include <iostream>

using namespace std;

// NOTE: Cycles start at 0.
void intensity2removezero(char *input_intensityfile,char *input_noisefile,char *output_intensityfile,char *output_noisefile) {

  // Intensity file reader
  IntensityReader input_intensity_file(input_intensityfile);
  IntensityReader input_noise_file    (input_noisefile);
  IntensityWriter output_intensity_file(output_intensityfile);
  IntensityWriter output_noise_file    (output_noisefile);

  for(;!input_intensity_file.eof();) {
   
    // An intensitysequence is a single line of the intensity file
    IntensitySequence is = input_intensity_file.get_next();
    IntensitySequence in = input_noise_file.get_next();

    if(!input_intensity_file.eof()) {
      bool badread=false;

      for(IntensitySequence::intensities_type::iterator i = is.intensities.begin();i != is.intensities.end();i++) {
        if((*i).anyoffedge()) badread=true;
      }
      
      if(!badread) {
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
