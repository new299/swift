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

template <class _prec>
inline std::ostream& operator<<(std::ostream& out, const Cluster<_prec> &rhs) {

  const vector<string> signal_ids = rhs.get_signal_ids();

  // Print intensities
  for(vector<string>::const_iterator i = signal_ids.begin();i != signal_ids.end();i++) {
    out << ">SIGNAL " << (*i) << " Position: " << rhs.get_position() << std::endl;
    int num_cycles = rhs.const_signal((*i)).size();
    for(int cycle=0;cycle < num_cycles;cycle++) {
      for(ReadIntensity<>::base_type base=0;base < ReadIntensity<_prec>::base_count;base++) {
        out << rhs.const_signal((*i))[cycle].get_base(base) << " ";
      }
      out << " ";
    }
    out << std::endl; 

    out << ">NOISE ESTIMATE : " << (*i) << " Position: " << rhs.get_position() << std::endl;
    num_cycles = rhs.const_noise((*i)).size();
    for(int cycle=0;cycle < num_cycles;cycle++) {
      for(ReadIntensity<>::base_type base=0;base < ReadIntensity<_prec>::base_count;base++) {
        out << rhs.const_noise((*i))[cycle].get_base(base) << " ";
      }
      out << " ";
    }
    out << std::endl;
  }
  
  
  const vector<string> sequence_ids = rhs.get_sequence_ids();
  
  // Print base calls
  for(vector<string>::const_iterator i = sequence_ids.begin();i != sequence_ids.end();i++) {
    out << ">SEQUENCE : " << (*i) << " Position: " << rhs.get_position() << std::endl;
    out << rhs.const_sequence(*i); // prints it's own linefeed
  }

  return out;
}
