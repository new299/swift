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

template <class _baset,class _qualt>
inline std::ostream& operator<<(std::ostream& out, const ScoredSequence<_baset,_qualt> &rhs) {

  for(typename ScoredSequence<_baset,_qualt>::sequence_type::const_iterator i=rhs.const_sequence().begin();i != rhs.const_sequence().end();i++) {
    out << ScoredSequence<_baset,_qualt>::base_name[(*i)];
  }

  out << std::endl;

  vector<int> q = rhs.phred_quality();
  for(typename vector<int>::const_iterator i=q.begin();i != q.end();i++) {
    out << (*i) << " "; // rescaled so highest quality is Q40, lowest is 1 in 4.
  }

  out << std::endl;

  return out;
}
