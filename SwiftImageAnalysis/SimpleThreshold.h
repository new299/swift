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

#ifndef SWIFTIMAGEANALYSIS_SIMPLETHRESHOLD
#define SWIFTIMAGEANALYSIS_SIMPLETHRESHOLD

#include <iostream>
#include <vector>
#include <math.h>
#include <stdexcept>
#include "SwiftImage.h"
#include "SwiftWindow.h"
#include "Timetagger.h"

using namespace std;

template <class _prec=uint16>
class SimpleThreshold {
public:

  SimpleThreshold(_prec threshold_in = 100) : threshold(threshold_in) {
  }

  SwiftImage<_prec> process(const SwiftImage<_prec> &source) {
  
    SwiftImage<_prec> dest = source;

    for(int n=0;n<source.image.size();n++) {
      if(dest.image[n] > threshold) dest.image[n] =1;
                               else dest.image[n] =0;
    }

    return dest;
  }

private:
  _prec threshold;         
};

#endif
