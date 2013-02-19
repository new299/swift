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

#ifndef SWIFTIMAGEANALYSIS_SWIFTIMAGEPOSITION
#define SWIFTIMAGEANALYSIS_SWIFTIMAGEPOSITION

/// This class representation a position on a SwiftImage
template <class _prec=int>
class SwiftImagePosition {

public:

  SwiftImagePosition() : x(0), y(0) {
  }

  SwiftImagePosition(_prec x_in,_prec y_in) : x(x_in), y(y_in) {
  }


  inline SwiftImagePosition<_prec> & operator=(const SwiftImagePosition &rhs) {
    x = rhs.x;
    y = rhs.y;
    return *this;
  }
  
  /// add positions
  inline SwiftImagePosition<_prec> operator+(const SwiftImagePosition &rhs) const {
    return SwiftImagePosition<_prec>(x+rhs.x,y+rhs.y);
  }
  
  /// subtract positions
  inline SwiftImagePosition<_prec> operator-(const SwiftImagePosition &rhs) const {
    return SwiftImagePosition<_prec>(x-rhs.x,y-rhs.y);
  }
  
  /// add positions
  inline SwiftImagePosition<_prec> operator+=(const SwiftImagePosition &rhs) {
    x=x+rhs.x;
    y=y+rhs.y;
  
    return (*this);
  }

  /// Are two positions the same
  inline bool operator==(SwiftImagePosition<_prec> rhs) {
    if((x == rhs.x) && (y == rhs.y)) return true;
                                else return false;
  }
  
  /// Are two positions the same
  inline bool operator!=(SwiftImagePosition<_prec> rhs) {
    return !((*this) == rhs);
  }

  //TODO: Add, and use, getters and setters
  _prec x; ///< X Position
  _prec y; ///< Y Position
};

#endif
