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

#ifndef SWIFTIMAGEANALYSIS_LOCALMAXIMA
#define SWIFTIMAGEANALYSIS_LOCALMAXIMA

#include <iostream>
#include <vector>
#include "SwiftImage.h"
#include <math.h>

using namespace std;


/// DEPRICATED THIS DOESN'T DO WHAT IT SHOULD
template <class _prec=uint16>
class LocalMaxima {
public:

  LocalMaxima() {
  }

  SwiftImage<_prec> process(const SwiftImage<_prec> &source) {
  
    SwiftImage<_prec> dest(source.image_width(),source.image_height());

    // For each pixel, examine a window_size by window_size window around the pixel,
    // Set the next pixel value to be 1 greater than the smallest value found.

    // Perform this left to right, top to bottom
    for(int x=0;x<source.image_width();x++) {
      for(int y=0;y<source.image_height();y++) {
        
        bool first=true;
        _prec max=0;
        for(int nx=-1;nx<=1;nx++) {
          for(int ny=-1;ny<=1;ny++) {

            // Not off the image, not this pixel
            if(((x+nx) >= 0) && ((x+nx) < source.image_width()) && (((y+ny) >= 0) && ((y+ny) < source.image_height())) && ((nx!=0) || (ny!=0))) {
              if((source(x+nx,y+ny) > max) || first) {
                max = source(x+nx,y+ny);
                first=false;
              }
            }
          }
        }
        
        // If a EDM pixel (non-zero) and greater than or equal to max, set it.
        if(source(x,y) != 0) {
          if(source(x,y) > max) dest(x,y) = 65534;
        }
      }
    }
    
    return dest;
  }

private:
  int window_size;
  vector<vector<bool> > mask;
};

#endif
