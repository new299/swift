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

#ifndef SWIFT_READ_IMAGELIST_H
#define SWIFT_READ_IMAGELIST_H

#include <iostream>
#include <fstream>

using namespace std;

vector<string> read_imagelist(char* imagelist_file) {
  vector<string> imagelist;

  ifstream input_file(imagelist_file);

  for(;!input_file.eof();) {
    string str;
    getline(input_file, str);
  
    if(!input_file.eof()) {
      imagelist.push_back(str);
    }
  }

  return str;
}

#endif
