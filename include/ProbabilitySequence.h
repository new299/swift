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


#ifndef SWIFT_PROBABILITYSEQUENCE_H
#define SWIFT_PROBABILITYSEQUENCE_H

#include <string>
#include <vector>
#include <math.h>
#include "BaseProbability.h"


using namespace std;

/// Represents a sequence with associated quality score
template<class _base_type=size_t,class _probability_prec=float>
class ProbabilitySequence {
public:

  typedef _base_type                        base_type;
  typedef _probability_prec                 probability_type;
  typedef string                            id_type;
  
  constexpr static const base_type base_a       = 0;    ///< Constant for base A
  constexpr static const base_type base_c       = 1;    ///< Constant for base C
  constexpr static const base_type base_g       = 2;    ///< Constant for base G
  constexpr static const base_type base_t       = 3;    ///< Constant for base T
  constexpr const static size_t    base_count   = 4;    ///< Number of bases
  
  typedef vector<BaseProbability<probability_type> > probability_sequence_type;

  static const std::string base_name[];  ///< string names for bases
  static const std::string base_descriptions[];  ///< string names for bases

  ProbabilitySequence() {
  }

  probability_sequence_type &sequence() {
    return m_sequence;
  }

  const probability_sequence_type &const_sequence() const {
    return m_sequence;
  }

  static string get_fast4_header() {
    ostringstream s;

    // Base format
    for(size_t n=0;n<base_count;n++) {
      s << "\\Base=" << n << "," << base_name[n] << "," << base_descriptions[n] << endl;
    }

    return s.str();
  }

  string get_fast4_string_sequence() const {
    
    ostringstream s;

    s << "@" << m_id << endl;

    for(size_t n=0;n<base_count;n++) {
      for(typename probability_sequence_type::const_iterator i=const_sequence().begin();i != const_sequence().end();i++) {
        s << setw(4) << (*i)[n] << " ";
      }
      s << endl;
    }

    return s.str();
  }

  string get_sequence_string() const {
    
    string s;
    for(typename probability_sequence_type::const_iterator i=const_sequence().begin();i != const_sequence().end();i++) {

      _probability_prec max=0;
      size_t max_idx       =0;
      for(size_t n=0;n<base_count;n++) {
        if((*i)[n] > max) {
          max = (*i)[n];
          max_idx = n;
        }
      }

      // don't call weird states
      //if((max_idx != base_complex) && (max_idx != base_end)) {
      s += base_name[max_idx];
      //}
    }

    return s;
  }

  const vector<int> phred_quality() const {
    vector<int> q;
    q.clear();

    for(typename probability_sequence_type::const_iterator i=const_sequence().begin();i != const_sequence().end();i++) {
      _probability_prec maxprb=0;
      for(size_t n=0;n<base_count;n++) {
        if((*i)[n] > maxprb) maxprb=(*i)[n];
      }

      q.push_back(static_cast<int>((30*maxprb)+1+0.5));
    }

    return q;
  }

  ProbabilitySequence<_base_type,_probability_prec> trim(int start,int end) const {
    ProbabilitySequence<_base_type,_probability_prec> trimmed;

    if(end == -1) end=m_sequence.size()-1;
    for(int i=start;i<=end;i++) {
      trimmed.m_sequence.push_back(m_sequence[i]);
    }

    trimmed.set_id(get_id());

    return trimmed;
  }



  id_type get_id() const {
    return m_id;
  }

  void set_id(const id_type &i) {
    m_id = i;
  }

  size_t size() const {
    return m_sequence.size();

  }

private:
  id_type                   m_id;       ///< Sequence ID
  probability_sequence_type m_sequence; ///< The sequence
};

template<class base_type,class qualityscore_type>
const std::string ProbabilitySequence<base_type,qualityscore_type>::base_name[] = {"A","C","G","T","INVALID"};

template<class base_type,class qualityscore_type>
const std::string ProbabilitySequence<base_type,qualityscore_type>::base_descriptions[] = {"Adenine","Cytosine","Guanine","Thymine","INVALID"};

#include "ProbabilitySequence.cpp"

#endif
