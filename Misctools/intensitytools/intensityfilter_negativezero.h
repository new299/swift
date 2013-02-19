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

#ifndef INTENSITYTOOLS_INTENSITYFILTER_NEGATIVEZERO_H
#define INTENSITYTOOLS_INTENSITYFITLER_NEGATIVEZERO_H

#include <vector>
#include <iostream>

using namespace std;

/// If anything is negative make it zero
class IntensityFilter_NegativeZero {
public:

  IntensityFilter_NegativeZero() {
  }

  IntensitySequence filter(IntensitySequence is) {
    // Walk across sequence anything that's negative becomes zero

    for(IntensitySequence::intensities_type::iterator i = is.intensities.begin();i != is.intensities.end();i++) {
      if( (*i).getbase(ReadIntensity::base_a) < 0) (*i).setbase(ReadIntensity::base_a,1);
      if( (*i).getbase(ReadIntensity::base_t) < 0) (*i).setbase(ReadIntensity::base_t,1);
      if( (*i).getbase(ReadIntensity::base_g) < 0) (*i).setbase(ReadIntensity::base_g,1);
      if( (*i).getbase(ReadIntensity::base_c) < 0) (*i).setbase(ReadIntensity::base_c,1);
    }

    return is;
  }
};

#endif
