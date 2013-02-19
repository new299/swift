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

#ifndef INTENSITYTOOLS_INTENSITYFILTER_PURITY_H
#define INTENSITYTOOLS_INTENSITYFITLER_PURITY_H

#include <vector>
#include <iostream>
#include "intensity.h"

using namespace std;

/// Filter based on purity
class IntensityFilter_Purity {
public:

  IntensityFilter_Purity(int num_bases_in=12,
                         double purity_threshold_in=0.6) 
                       : num_bases(num_bases_in),
                         purity_threshold(purity_threshold_in) {
  }

  IntensitySequence filter(IntensitySequence is) {
    
    if(is.intensities.size() > (num_bases+1))
    if(is.min_purity(0,num_bases) > purity_threshold) {
      is.valid = true;
    } else is.valid = false;
    
    return is;
  }

  int num_bases;
  double purity_threshold;
};

#endif
