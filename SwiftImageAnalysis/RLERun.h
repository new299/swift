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

#ifndef SWIFTIMAGEANALYSIS_RLERUN
#define SWIFTIMAGEANALYSIS_RLERUN

/// This class represents a Run Length Encoded region of an image.
template <class _prec=int>
class RLERun {
  public:

  RLERun(_prec x_in,_prec y_in,_prec length_in) : pos(x_in,y_in),
                                                  length(length_in) {
  }


  template<class _iprec>
  _iprec max_pixel(const SwiftImage<_iprec> &i,bool &onimage) const {
    _prec y=pos.y;
    bool first=true;
    _iprec maxval=0;
    onimage=false;
    for(_prec x=pos.x;x < pos.x+length;x++) {
      if(i.onimage(x,y)) {
        _iprec val = i(x,y);
        if((val > maxval) || first) {
          maxval = val;
          onimage=true;
          first=false;
        }
      }
    }

    return maxval;
  }



  SwiftImagePosition<_prec> pos;
  _prec length;
};

#endif
