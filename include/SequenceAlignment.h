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

#ifndef SWIFT_SEQUENCEALIGNMENT_H
#define SWIFT_SEQUENCEALIGNMENT_H

template<class _prec=int>
class SequenceAlignment {
public:

  SequenceAlignment() {
  }

  SequenceAlignment(_prec contig_in,
                    _prec position_in,
                    _prec score_in,
                    vector<bool> matchstring_in,
                    bool unique_in) :
                    contig(contig_in),
                    position(position_in),
                    score(score_in),
                    matchstring(matchstring_in),
                    unique(unique_in) {
  }
  
  _prec contig;
  _prec position;
  _prec score;
  vector<bool> matchstring;
  bool  unique;
};

#endif
