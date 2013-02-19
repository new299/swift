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

#include "SwiftImage.h"

using namespace std;

int main(int argc,char **argv) {

  if(argc < 2) {
    cout << "imghist <input tiff>" << endl;
    return 0; 
  }

  SwiftImage<> s(argv[1]);

  vector<int> hist(70000,0);

  for(size_t x=s.min_x();x<s.max_x();x++) {
    for(size_t y=s.min_y();y<s.max_y();y++) {
      int val = s(x,y);
      hist[val]++;
    }
  }

  for(int n=0;n<hist.size();n++) {
    cout << n << " " << hist[n] << endl;
  }

  return 0;
}
