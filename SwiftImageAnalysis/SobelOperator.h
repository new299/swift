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

#ifndef SWIFTIMAGEANALYSIS_SOBELOPERATOR
#define SWIFTIMAGEANALYSIS_SOBELOPERATOR

#include <iostream>
#include <vector>
#include "SwiftImage.h"
#include <math.h>

using namespace std;

/// Class to apply a sobel operator to an image. Works, but not currently used.
template <class _prec=uint16>
class SobelOperator {
public:

  SobelOperator() {
  }

  // Applies Sobel Kernel to an image
  SwiftImage<_prec> process(const SwiftImage<_prec> &source) {
    SwiftImage<_prec> dest(source.image_width(),source.image_height());
    
    for(int x=1;x < (source.image_width()-1);x++) {
      for(int y=1;y < (source.image_height()-1);y++) {
        _prec p1 = source(x-1,y-1);
        _prec p2 = source(x  ,y-1);
        _prec p3 = source(x+1,y-1);
        _prec p4 = source(x-1,y  );
        // _prec p5 = source(x  ,y  );
        _prec p6 = source(x+1,y  );
        _prec p7 = source(x-1,y+1);
        _prec p8 = source(x  ,y+1);
        _prec p9 = source(x+1,y+1);

        int i = (p1+2*p2+p3)-(p7+2*p8+p9);
        int j = (p3+2*p6+p9)-(p1+2*p4+p7);

        if(i < 0) i = 0-i;
        if(j < 0) j = 0-j;

        dest(x,y) = i+j;
      }
    }

    return dest;
  }

private:
};

#endif
