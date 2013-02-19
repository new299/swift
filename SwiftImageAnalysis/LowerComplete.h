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

#ifndef SWIFTIMAGEANALYSIS_LOWERCOMPLETE
#define SWIFTIMAGEANALYSIS_LOWERCOMPLETE

#include <iostream>
#include <vector>
#include <deque>
#include "SwiftImage.h"
#include "SwiftImagePosition.h"
#include <math.h>

using namespace std;

template <class _prec=uint16>
class LowerComplete {
public:

  LowerComplete() {
  }

  bool neighbour_lower(const SwiftImage<_prec> &source,int x,int y) {
    // for all neighbours
    for(int nx=-1;nx<=1;nx++) {
      for(int ny=-1;ny<=1;ny++) {
        if(!((nx == 0) && (ny == 0))) {
          
          if(source.onimage(x,y) && source.onimage(x+nx,y+ny)) {
            if(source(x+nx,y+ny) < source(x,y)) return true; //< changed to <=
          }
        }
      }
    }
    
    return false;
  }
  
  SwiftImage<_prec> process(const SwiftImage<_prec> &source) {
 
    // needs to be signed
    SwiftImage<int16> outimg(0,0);
    outimg = source;
    SwiftImage<_prec> outimg_real = source;

    deque<SwiftImagePosition<int> > queue;

    SwiftImagePosition<int> end_marker(-1,-1);

    // populate queue
    // cout << "queue: ";
    for(int x=source.min_x();x<source.max_x();x++) {
      for(int y=source.min_y();y<source.max_y();y++) {
        outimg(x,y) = 0;
        if(neighbour_lower(source,x,y)) {
          outimg(x,y) = -1;
          queue.push_back(SwiftImagePosition<int>(x,y));
          // cout << x << "," << y << " ";
        }
      }
    }
    // cout << endl;

    queue.push_back(end_marker);
    int dist=1;

    while(queue.size() > 0) {
      SwiftImagePosition<int> c = static_cast<SwiftImagePosition<int> >(queue.front());
      queue.pop_front();

      if(c == end_marker) {
        if(queue.size() > 0) {
          if(queue[queue.size()-1] != end_marker) {
            queue.push_back(end_marker);
            dist++;
          }
          // cout << "intermed:" << endl;
          // outimg.dump(cout);
        }
      } else {
        outimg(c) = dist;
        // cout << "outimg: " << c.x << "," << c.y << "=" << dist << endl;
        
        // for all neighbours
        for(int nx=-1;nx<=1;nx++) {
          for(int ny=-1;ny<=1;ny++) {
            if(!((nx == 0) && (ny == 0))) {
              
              // If they are the same height add neighbours to the queue
              if(source.onimage(c.x+nx,c.y+ny) && source.onimage(c.x,c.y) && outimg.onimage(c.x+nx,c.y+ny)) {
                if((source(c.x+nx,c.y+ny) == source(c.x,c.y)) && (outimg(c.x+nx,c.y+ny) == 0)) {
                  queue.push_back(SwiftImagePosition<int>(c.x+nx,c.y+ny));
                  // cout << "pushing back: " << c.x+nx << "," << c.y+ny << endl;
                  outimg(c.x+nx,c.y+ny) = -1;
                }
              }
            }
          }
        }
      }
    }

    for(int x=outimg.min_x();x<outimg.max_x();x++) {
      for(int y=outimg.min_y();y<outimg.max_y();y++) {
        if(outimg(x,y) != 0) {
          outimg(x,y) = dist*source(x,y)+outimg(x,y)-1;
        } else {
          // outimg(x,y) = dist*source(x,y);
        }
      }
    }
    
    outimg_real = outimg;
    
    return outimg_real;
    
  }


};

#endif
