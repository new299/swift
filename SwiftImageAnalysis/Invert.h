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

#ifndef SWIFTIMAGEANALYSIS_INVERT
#define SWIFTIMAGEANALYSIS_INVERT

#include <iostream>
#include <vector>
#include "SwiftImage.h"
#include <math.h>

using namespace std;

/// Invert an image
template <class _prec=uint16>
class Invert {
public:

  Invert() {
  }

  SwiftImage<_prec> process(const SwiftImage<_prec> &source) {

    _prec maxval = source.max();

    SwiftImage<_prec> dest = source;

    // For each pixel, examine a window_size by window_size window around the pixel,
    // Set this pixel to the minimum value found.

    for(int x=source.min_x();x<source.max_x();x++) {
      for(int y=source.min_y();y<source.max_y();y++) {
        dest(x,y) = maxval-source(x,y); 
      }
    }
    
    return dest;
  }
};

#endif
