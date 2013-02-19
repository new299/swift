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

#ifndef SWIFT_CLUSTERPOSITION_H
#define SWIFT_CLUSTERPOSITION_H

#include <string>

template <class _prec=int> 
class ClusterPosition {
  public:

  _prec x;
  _prec y;

   ClusterPosition() {
   }
   
   ClusterPosition(_prec x_in,_prec y_in) : x(x_in), y(y_in) {
   }


  void set(_prec x_in,_prec y_in) {
    x=x_in;
    y=y_in;
  }

  std::string as_string() const {
    std::string s;
    s += stringify(x);
    s += " ";
    s += stringify(y);

    return s;
  }
};

#include "ClusterPosition.cpp"

#endif
