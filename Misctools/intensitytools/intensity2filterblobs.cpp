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

#include "intensity2filterblobs.h"
#include "stringify.h"

using namespace std;

int main(int argc,char **argv) {

  if(argc < 3) {
    cout << "intensity2filterblobs <input intensity file> <input noise file> <output intensity file> <output noise file> <no bad cycles allowed>" << endl;
    cout << "remove reads where intensity is constant or increases across cycle" << endl;
    return 0; 
  }

  intensity2filterblobs(argv[1],argv[2],argv[3],argv[4],convertTo<int>(argv[5]));

  return 0;
}
