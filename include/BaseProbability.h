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


#ifndef SWIFT_BASEPROBABILITY_H
#define SWIFT_BASEPROBABILITY_H

#include <string>
#include <vector>
#include <math.h>
#include <sstream>

using namespace std;

/// Represents a sequence with associated quality score
template<class _base_type=size_t,class _probability_prec=float>
class BaseProbability {
public:

  typedef _base_type                        base_type;
  typedef _probability_prec                 probability_type;
  
  constexpr static const base_type base_a       = 0;    ///< Constant for base A
  constexpr static const base_type base_c       = 1;    ///< Constant for base C
  constexpr static const base_type base_g       = 2;    ///< Constant for base G
  constexpr static const base_type base_t       = 3;    ///< Constant for base T
  constexpr const static size_t    base_count   = 4;    ///< Number of bases
  

  static const std::string base_name[];  ///< string names for bases
  static const std::string base_descriptions[];  ///< string names for bases

  BaseProbability() {
  }

  inline _probability_prec &operator[] (size_t idx) {
    return m_probability[idx];
  }
  
  inline const _probability_prec &operator[] (size_t idx) const {
    return m_probability[idx];
  }

  static string get_fast4_header() {
    ostringstream s;

    // Base format
    for(size_t n=0;n<base_count;n++) {
      s << "\\Base=" << n << "," << base_name[n] << "," << base_descriptions[n] << endl;
    }

    return s.str();
  }

private:
  _probability_prec m_probability[base_count]; ///< The probability
};

template<class base_type,class qualityscore_type>
const std::string BaseProbability<base_type,qualityscore_type>::base_name[] = {"A","C","G","T","INVALID"};

template<class base_type,class qualityscore_type>
const std::string BaseProbability<base_type,qualityscore_type>::base_descriptions[] = {"Adenine","Cytosine","Guanine","Thymine","INVALID"};

#endif
