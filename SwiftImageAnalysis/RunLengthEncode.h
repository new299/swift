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

#ifndef SWIFTIMAGEANALYSIS_RUNLENGTHENCODE
#define SWIFTIMAGEANALYSIS_RUNLENGTHENCODE

#include <iostream>
#include <vector>
#include "SwiftImage.h"
#include "RLERun.h"
#include <math.h>

using namespace std;

template <class _prec=uint16>
class RunLengthEncode {
public:
  
  RunLengthEncode() {
  }

  // Convert an image in to a set of RLERuns
  vector<RLERun<> > process(const SwiftImage<_prec> &source) {
   
    vector<RLERun<> > runs;

    bool in_run      = false;
    int  run_start_x = 0;
    int  run_length  = 0;

    for(int y=source.min_y();y<source.max_y();y++) {
      for(int x=source.min_x();x<source.max_x();x++) {
        if(source.onimage(x,y)) { 
          if(source(x,y) > 0) {
            if(in_run) {
              run_length++;
            } else {
              in_run      = true;
              run_start_x = x;
              run_length  = 1;
            }
          } else {
            if(in_run) {
              runs.push_back(RLERun<int>(run_start_x,y,run_length));
            }
            in_run=false;
          }
        } else {}//cerr << "not onimage: " << x << "," << y << endl;}
      }

      if(in_run) {
        runs.push_back(RLERun<int>(run_start_x,y,run_length));
      }
      in_run=false;
    }

    return runs;
  }

private:
};

#endif
