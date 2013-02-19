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
#include "ProbabilitySpecification.h"

using namespace std;

/// Represents a sequence with associated quality score

template<class _base_probability_type>
class ProbabilitySequence {
public:

  typedef string                         id_type;
  typedef vector<_base_probability_type> probability_sequence_type;

  ProbabilitySequence() {

  }

  void set_specification(string fast4header) {

  }

  probability_sequence_type &sequence() {
    return m_sequence;
  }

  const probability_sequence_type &const_sequence() const {
    return m_sequence;
  }

  string string_sequence() const {
    string s;

    for(int n=0;n<m_sequence.size();n++) {
      string sn = m_sequence[n].short_name();
      s += sn;
    }

    return s;
  }

  ProbabilitySequence<_base_probability_type> trim(int start,int end) const {
    ProbabilitySequence<_base_probability_type> trimmed;

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

#endif
