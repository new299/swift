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

#ifndef SWIFTIMAGEANALYSIS_DSETS
#define SWIFTIMAGEANALYSIS_DSETS

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

class DSets {
public:

   DSets(int dsets_size) {
     parents.clear();
     parents.insert(parents.begin(),dsets_size,-1);
   }

   void makeparent(int c,int l) {
     // cout << "makeparent called: " << c << "," << l << endl;
     parents[c] = l;
   }

   int getparent(int c) {
     
     if(parents[c] == -1) {
       return c;
     } else {
       int r = getparent(parents[c]);
       parents[c] = r;
       return parents[c];
     }
   
   }

   void dump(ostream &out) {
     for(unsigned int n=0;n<parents.size();n++) {
       out << n << ": " << parents[n] << endl;
     }
   }

private:
  vector<int> parents;
};

#endif
