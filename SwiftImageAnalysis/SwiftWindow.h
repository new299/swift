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

#ifndef SWIFTWINDOW_H_
#define SWIFTWINDOW_H_

#include <stdexcept>
#include <iomanip>
#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include "SwiftImage.h"

using namespace std;

template <class _prec=uint16>
class SwiftWindow {

public:

  typedef enum min_max_t {find_min, find_max};

  SwiftWindow (int window_size=6) : m_window_size(window_size) {};
  
  SwiftImage<_prec> window_square (const SwiftImage<_prec> &source, min_max_t which) {

    int x_dim = source.image_width();
    int y_dim = source.image_height();
    
    const vector<_prec> &from = source.get_image();

    vector<_prec> accum1 (x_dim*y_dim, 0);
    vector<_prec> accum2 (x_dim*y_dim, 0);

    // We're going to find the max and min values in a window around each element in two
    // passes: First in the x direction, finding the max/min of the adjacent row elements.
    // Then we process the result of that in the y direction, to find the maximum of the
    // maximums (or minimum of the minimums). window_pass returns its result in transposed
    // order to facilitate this.

    // We'll use the SwiftImage convention for indexing a vector like a 2-d array: The x
    // coordinate varies the fastest -- so two adjacent elements in the vector are adjacent
    // in a row of the conceptual 2-d array. For transposed results, the opposite is true.
    
    // Because 2-d array may not be square, we need to tell window_pass its x-dimension.
    // The y-dimension can then be inferred from the vector size.
    
    window_pass (from, accum1, x_dim, which);     // returns transposed result
    window_pass (accum1, accum2, y_dim, which);   // transposed again -- back to normal!
    
    SwiftImage<_prec> dest(x_dim, y_dim);
    vector<_prec> &img = dest.get_image();
    img.clear();
    
    // In a SwiftImage, the row (x) coordinate varies the fastest. 
    
    for (int y=0; y<y_dim; y++) {
      for (int x=0; x<x_dim; x++) {
        int index = y*x_dim+x;
        img.push_back (accum2[index]);
      }
    }

    return dest;
    
  }
    
  // Print vector as 2-d array. row_len is size of fastest-varying index.
  
  void print_vec (const vector<_prec> &vec, int row_len,
      int start_x=0, int start_y=0, int square=-1) {
     
    int from_x = (start_x >= 0) ? start_x : 0;
    int from_y = (start_y >= 0) ? start_y : 0;
    
    int x_dim = row_len;
    int y_dim = vec.size() / x_dim;
    
    int to_x = (square == -1) ? x_dim : from_x+square;
    int to_y = (square == -1) ? y_dim : from_y+square;
    
    printf ("x: %4i  y: %4i  size: %8ld\n", x_dim, y_dim, vec.size());
    
    for (int y=from_y; y<to_y; y++) {
      for (int x=from_x; x<to_x; x++) {
        int index = y*x_dim+x;  
        cout << right << setw(3) << vec[index] << " ";  
      }
      cout << endl;
    }
    cout << "\n" << endl;

  }

private:
  
  void window_pass (const vector<_prec> &from, vector<_prec> &to, int row_len, min_max_t which) {

    int x_dim = row_len;
    int y_dim = from.size() / x_dim;
    
    for (int y=0; y<y_dim; y++) {

      list<_prec> Q;
      Q.assign (m_window_size, from[y*x_dim]);    // preload trailing values with zeroth element of row

      int lead = 0;                               // x index of leading value
      for ( ; lead <= m_window_size; lead++) {    // load leading values
        int index = y*x_dim+lead;
        Q.push_back (from[index]);
      }

      _prec cur_best = (which == find_max) ? max_of(Q) : min_of(Q);

      // At this point, Q contains 2*window+1 values: window_size phoney trailing
      // elements (copies of element 0,y), the target element (0,y in this case),
      // plus the next window_size values from the array. At each loop iteration,
      // we set the current element to cur_best, which is the max (or min) of the
      // values in Q. Then we pop the oldest value off Q, and push a new one on.
      // If the new pushed value is >= cur_best, it becomes the new cur_best. Else
      // if the popped value was the maximum, we need a new maximum, which requires
      // examining all the values in Q. If neither is true, on to the next element,
      // leaving cur_best unchanged.

      for (int x=0; x<x_dim; x++, lead++) {  // process one row

        to[x*y_dim+y] = cur_best;            // store into array transposed !!!!

        _prec pop_off = Q.front();
        Q.pop_front();                       // doesn't return the popped item (aargh)

        int index = y*x_dim+lead;
        if (lead >= x_dim) {
          index = y*x_dim+x_dim-1;           // index of last element of row;
        }
        _prec push_on = from[index];
        Q.push_back (push_on);

        if (which == find_max) {
          
          if (push_on >= cur_best) {
            cur_best = push_on;                   // we've found a new max
          } else if (pop_off >= cur_best) {       // old max is gone, need a new one
            cur_best = max_of(Q);
          }                                      // else old cur_best is still good
          
        } else {
          
          if (push_on <= cur_best) {
            cur_best = push_on;
          } else if (pop_off >= cur_best) {
            cur_best = min_of(Q);
          }
          
        }
        
      }

    }
    
//////    print_vec (to, y_dim);                     // remember, to is transposed

  }
  
  _prec max_of (const list<_prec> &Q) {

    typename list<_prec>::const_iterator iter;

    _prec cur_max = Q.front();
    for (iter=Q.begin(); iter != Q.end(); iter++) {
      if (*iter > cur_max) {
        cur_max = *iter;
      }
    }
    
    return cur_max;
    
  }
  
  _prec min_of (const list<_prec> &Q) {

    typename list<_prec>::const_iterator iter;

    _prec cur_min = Q.front();
    for (iter=Q.begin(); iter != Q.end(); iter++) {
      if (*iter < cur_min) {
        cur_min = *iter;
      }
    }
    
    return cur_min;
    
  }
  
  void print_Q (const list<_prec> &Q) {

    typename list<_prec>::const_iterator iter;

    for (iter=Q.begin(); iter != Q.end(); iter++) {
      cout << right << setw(3) << *iter << " ";
    }
    cout << endl;
    
  }
  
  int m_window_size;

};
  
#endif /*SWIFTWINDOW_H_*/
