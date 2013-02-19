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

#ifndef SWIFTIMAGEANALYSIS_WATERSHED
#define SWIFTIMAGEANALYSIS_WATERSHED

#include <iostream>
#include <vector>
#include "SwiftImage.h"
#include <math.h>
#include "LowerComplete.h"
#include "DSetsWatershed.h"

using namespace std;


/// This class performs Watershed segmentation using the Disjoint set algorithm.
template <class _prec=uint16>
class Watershed {
public:

  Watershed() {
  }

  bool check_adj(const SwiftImage<_prec> &i,vector<int> &l,_prec min_neighbour,DSetsWatershed &s,int x,int y,int x_off,int y_off) {
    if(i.onimage(x+x_off,y+y_off)) {
      if(i(x+x_off,y+y_off) == min_neighbour) {
        if(s.getparents(i.get_index(x+x_off,y+y_off)).size() != 0) {
          if(s.getparents(static_cast<int>(i.get_index(x+x_off,y+y_off)))[0] != static_cast<int>(i.get_index(x,y))) {
            l.push_back(i.get_index(x+x_off,y+y_off));
            return true;
          } else return false;
        } else {
          l.push_back(i.get_index(x+x_off,y+y_off));
          return true;
        }
      } else return false;
    
    } else return false;

    return false;
  }

  vector<int> lowest_neighbours(const SwiftImage<_prec> &i,DSetsWatershed &s,int x,int y) {

    bool first=true;
    int min_neighbour=0;

    for(int nx=-1;nx<=1;nx++) {
      for(int ny=-1;ny<=1;ny++) {
        
        if(i.onimage(x+nx,y+ny)) {

          if((i(x+nx,y+ny) < min_neighbour) || first) {
            first=false;
            min_neighbour=i(x+nx,y+ny);
            // cout << "x+nx: " << x+nx << endl;
            // cout << "y+ny: " << y+ny << endl;
            // cout << "nx: " << nx << endl;
            // cout << "ny: " << ny << endl;
          }
        }
      }
    }
    
    vector<int> l;
    if(i.onimage(x,y)) {
      if(min_neighbour < i(x,y)) {
        for(int nx=-1;nx<=1;nx++) {
          for(int ny=-1;ny<=1;ny++) {
            if(i.onimage(x+nx,y+ny)) {
              
              if(i(x+nx,y+ny) == min_neighbour) {
                // cout << "LOWEST NEIGHBOUR: " << x+nx << "," << y+ny << endl;
                l.push_back(i.get_index(x+nx,y+ny));
              }
            }
          }
        }
      } else {
        l.clear();

        // If on a plateau move anticlockwise if possible (to find canonical pixel).
        // Don't return from where we came from.
        
        if(min_neighbour == i(x,y)) {
          // This is why we can't have nice things... please rewrite this code.
          bool add=false;

          int nx,ny;
          // cout << "what what" << endl;
          if(add==false) {nx=0 ; ny=-1 ; add=check_adj(i,l,min_neighbour,s,x,y,nx,ny);}
          if(add==false) {nx=-1; ny=-1 ; add=check_adj(i,l,min_neighbour,s,x,y,nx,ny);}
          if(add==false) {nx=-1; ny=0 ;  add=check_adj(i,l,min_neighbour,s,x,y,nx,ny);}
          if(add==false) {nx=-1; ny=1;   add=check_adj(i,l,min_neighbour,s,x,y,nx,ny);}
          if(add==false) {nx=0 ; ny=1;   add=check_adj(i,l,min_neighbour,s,x,y,nx,ny);}
          if(add==false) {nx=1 ; ny=1;   add=check_adj(i,l,min_neighbour,s,x,y,nx,ny);}
          if(add==false) {nx=1 ; ny=0 ;  add=check_adj(i,l,min_neighbour,s,x,y,nx,ny);}
          if(add==false) {nx=1 ; ny=-1 ; add=check_adj(i,l,min_neighbour,s,x,y,nx,ny);}
          // cout << "what what end" << endl;

        } 
      }
    }

    return l;
  }

  SwiftImage<_prec> process(const SwiftImage<_prec> &source) {
  
    SwiftImage<_prec> dest = source;
    LowerComplete<_prec> lcomp;
    dest = lcomp.process(source);

    image_width=dest.image_width();
    image_height=dest.image_height();

    DSetsWatershed sets(source.max_index());

    //1. Make EDM image lower complete

    // For all pixels, build the Watershed forest
    for(int x=source.min_x();x<source.max_x();x++) {
      for(int y=source.min_y();y<source.max_y();y++) {

        vector<int> l = lowest_neighbours(dest,sets,x,y);
        if(l.size() == 0) {
          sets.add(source.get_index(x,y)); // add if it doesn't exist
        } else {
          for(vector<int>::iterator i=l.begin();i != l.end();i++) {
            //make c a child of lowestneighbour
            sets.makeparent(source.get_index(x,y),(*i));
          }
        }
      }
    }

    // Detect watershed pixels
    for(int x=source.min_x();x<source.max_x();x++) {
      for(int y=source.min_y();y<source.max_y();y++) {
        vector<int> p = sets.getparents(source.get_index(x,y));

        if(p.size() != 1) dest(x,y) = 0; else dest(x,y) = 1;
      }
    }

    return dest;
  }

private:
  int image_width;      ///< Width of the image being processed
  int image_height;     ///< Height of the image being processed
};

#endif
