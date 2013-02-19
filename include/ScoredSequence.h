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

#ifndef SWIFT_SCOREDSEQUENCE_H
#define SWIFT_SCOREDSEQUENCE_H

#include <string>
#include <vector>
#include <math.h>

using namespace std;

/// Represents a sequence with associated quality score
template<class _base_type=int,class _qualityscore_type=double>
class ScoredSequence {
public:

  typedef _base_type                base_type;
  typedef _qualityscore_type        qualityscore_type;
  typedef string                    id_type;
  typedef vector<base_type>         sequence_type;
  typedef vector<qualityscore_type> quality_type;
  
  static const base_type base_a = 0;          ///< Constant for base A
  static const base_type base_c = 1;          ///< Constant for base C
  static const base_type base_g = 2;          ///< Constant for base G
  static const base_type base_t = 3;          ///< Constant for base T
  static const base_type base_n = 4;          ///< Constant for base N/Unknown
  const static base_type base_count = 5;      ///< Number of bases

  static const std::string base_name[];       ///< string names for bases

  ScoredSequence() {
  }

  void set_quality(int new_quality) {
    quality.clear();
    quality.insert(quality.begin(),sequence.length(),new_quality);
  }
  
  sequence_type &sequence() {
    return m_sequence;
  }

  quality_type &quality() {
    return m_quality;
  }
  
  const sequence_type &const_sequence() const {
    return m_sequence;
  }

  const quality_type &const_quality() const {
    return m_quality;
  }

  string get_sequence_string() const {
    
    string s;
    for(typename sequence_type::const_iterator i=const_sequence().begin();i != const_sequence().end();i++) {
      s += base_name[(*i)];
    }

    return s;
  }

  id_type get_id() const {
    return id;
  }

  bool set_id(const id_type &i) {
    id = i;
    return true;
  }

  const vector<int> phred_quality() const {
    vector<int> q;
    q.clear();

    for(typename quality_type::const_iterator i=m_quality.begin();i != m_quality.end();i++) {
      q.push_back(static_cast<int>((34*(*i))+6+0.5));
    }

    return q;
  }
  
  ScoredSequence<_base_type,_qualityscore_type> trim(int start,int end) {
    ScoredSequence<_base_type,_qualityscore_type> trimmed;
   
    if(end == -1) end=m_sequence.size()-1; 
    for(int i=start;i<=end;i++) {
      trimmed.m_sequence.push_back(m_sequence[i]);
      trimmed.m_quality.push_back(m_quality[i]);
    }

    trimmed.id = id;
    
    return trimmed;
  }

private:
  id_type       id;          ///< Sequence ID
  sequence_type m_sequence;  ///< The sequence
  quality_type  m_quality;   ///< The associated quality score
};

template<class base_type,class qualityscore_type>
const std::string ScoredSequence<base_type,qualityscore_type>::base_name[] = {"A","C","G","T","N","INVALID"};

#include "ScoredSequence.cpp"

#endif
