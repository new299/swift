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

#ifndef SIMPLECALLER_INTENSITYFILE
#define SIMPLECALLER_INTENSITYFILE

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "intensity.h"

using namespace std;

class IntensityWriter {
  
  string   output_filename;
  ofstream output_file;

public:

  IntensityWriter() {
  }

  IntensityWriter(string filename) : output_filename(filename) {
    output_file.open(filename.c_str());
  }

  void open(string filename) {
    output_filename = filename;
    output_file.open(filename.c_str());
  }

  void write(const IntensitySequence &is) {
    output_file << is.lane << "\t" << is.tile << "\t" << is.x << "\t" << is.y;
    for(IntensitySequence::intensities_type::const_iterator i = is.intensities.begin();i != is.intensities.end();i++) {
      output_file << "\t";
      output_file << (*i).getbase(ReadIntensity::base_a) << " ";
      output_file << (*i).getbase(ReadIntensity::base_c) << " ";
      output_file << (*i).getbase(ReadIntensity::base_g) << " ";
      output_file << (*i).getbase(ReadIntensity::base_t);
    }
    output_file << endl;
  }

  void close() {
    output_file.close();
  }
};

#endif
