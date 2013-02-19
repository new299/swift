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

#ifndef INTENSITYTOOLS_INTENSITYFILTER_CONSTANT_H
#define INTENSITYTOOLS_INTENSITYFITLER_CONSTANT_H

#include <vector>
#include <iostream>
#include "intensity.h"

using namespace std;

class IntensityFilter_Constant {
public:

  IntensityFilter_Constant(double const_a_in,
                           double const_c_in,
                           double const_g_in,
                           double const_t_in) :
                           const_a(const_a_in),
                           const_c(const_c_in),
                           const_g(const_g_in),
                           const_t(const_t_in) {
  }

  IntensitySequence filter(IntensitySequence is) {
    // Walk across sequence anything that's negative becomes zero

    for(IntensitySequence::intensities_type::iterator i = is.intensities.begin();i != is.intensities.end();i++) {
      (*i).setbase(ReadIntensity::base_a,(*i).getbase(ReadIntensity::base_a)+const_a);  
      (*i).setbase(ReadIntensity::base_c,(*i).getbase(ReadIntensity::base_c)+const_c);  
      (*i).setbase(ReadIntensity::base_g,(*i).getbase(ReadIntensity::base_g)+const_g);  
      (*i).setbase(ReadIntensity::base_t,(*i).getbase(ReadIntensity::base_t)+const_t);  
    }

    return is;
  }

  double const_a;
  double const_c;
  double const_g;
  double const_t;
};

#endif
