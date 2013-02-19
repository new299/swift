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
#include "ProbabilitySpecification.h"

using namespace std;

/// Represents a sequence with associated quality score
template<int   _probability_id,
         int   _base_count      =4,
         class _probability_prec=float
        >
class BaseProbability {
public:

  typedef _probability_prec         probability_type;
  static ProbabilitySpecification m_probability_spec;
  static const int base_count = _base_count;

  BaseProbability(int c) {
  }
  
  static ProbabilitySpecification &get_probability_specification() {
    
    return m_probability_spec; // *(new ProbabilitySpecification()); // m_probability_spec;
  }

  inline _probability_prec &operator[] (size_t idx) {
    return m_probability[idx];
  }
  
  inline const _probability_prec &operator[] (size_t idx) const {
    return m_probability[idx];
  }

  inline const string& short_name() const {
    _probability_prec maxval = m_probability[0];
    size_t maxidx = 0;

    for(int n=0;n<_base_count;n++) {
      if(m_probability[n] >= maxval) {
        maxval = m_probability[n];
        maxidx = n;
      }
    }
    
    return m_probability_spec.get_name(maxidx);
  }

  void set_base(size_t idx,const _probability_prec val) {
    m_probability[idx] = val;
  }

  _probability_prec max_probability() {
    _probability_prec maxval = m_probability[0];

    for(int n=0;n<_base_count;n++) {
      if(m_probability[n] >= maxval) {
        maxval = m_probability[n];
      }
    }
   
    return maxval; 
  }

private:
  _probability_prec m_probability[_base_count]; ///< The probability
};

template<int _probability_id,int _base_count,class _probability_prec>
ProbabilitySpecification BaseProbability<_probability_id,_base_count,_probability_prec>::m_probability_spec;

#endif
