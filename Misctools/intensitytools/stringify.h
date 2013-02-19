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

#ifndef STRINGIFY_H
#define STRINGIFY_H

#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>

class BadConversion : public std::runtime_error {
public:
 BadConversion(const std::string& s)
      : std::runtime_error(s) {}
};

template<class _type>
inline std::string stringify(_type x)
{
  std::ostringstream o;
  if (!(o << std::fixed << x))
    throw BadConversion("stringify()");
  return o.str();
}

template<class _type>
inline _type convertTo(const std::string& s)
{
  std::istringstream i(s);
  _type x;
  if (!(i >> x))
    throw BadConversion("convertTo(\"" + s + "\")");
  return x;
}

#endif

