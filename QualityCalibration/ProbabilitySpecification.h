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


#ifndef SWIFT_PROBABILITYSPECIFICATION_H
#define SWIFT_PROBABILITYSPECIFICATION_H

#include <string>
#include <vector>
#include <map>
#include <math.h>

using namespace std;

/// Represents a sequence with associated quality score
class ProbabilitySpecification {
public:

  ProbabilitySpecification() {
  }

  size_t add_base(std::string new_name,
                  std::string new_description) {
    base_names.push_back(new_name);
    base_descriptions.push_back(new_description);
    base_index[new_name] = base_names.size()-1;
    
    return base_names.size()-1;
  }

  void load_fast4_header(std::string header) {
  }

  void add_specification(size_t index,
                         string shortname,
                         string longname) {
  
    base_names.push_back(shortname);
    base_descriptions.push_back(longname);
    base_index[shortname] = index;
  }
  
  // FIXME: only return value if contained in map, else return something else/throw exception.
  //        modify so does not change map
  size_t get_index(std::string lookup_name) {
    return base_index[lookup_name];
  }

  const std::string &get_description(size_t lookup_idx) const {
    return base_descriptions[lookup_idx];
  }

  // FIXME: modify so const
  const std::string &get_description(std::string lookup_name) {
    return base_descriptions[get_index(lookup_name)];
  }

  const std::string &get_name(size_t idx) const {
    return base_names[idx];
  }

  size_t get_probability_count() {
    return base_names.size();
  }

  void clear() {
    base_names.clear();
    base_index.clear();
    base_descriptions.clear();
  }

private:
  vector<std::string>             base_names;         ///< string names for bases
  map<std::string,size_t>         base_index;         ///< lookup index for names (inverse on base_name)
  vector<std::string>             base_descriptions;  ///< string names for bases
};

#endif
