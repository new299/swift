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

#include "intensityfilter_constant.h"
#include "stringify.h"
#include "intensity.h"
#include "intensityreader.h"
#include "intensitywriter.h"

using namespace std;

int main(int argc,char **argv) {

  if(argc < 6) {
    cout << "intensity2addconst <input intensity file> <output intensity file> <a const> <c const> <g const> <t const>" << endl;
    return 0; 
  }

  const char* input_intensity_filename  = argv[1];
  const char* output_intensity_filename = argv[2];

  double a_const = convertTo<double>(argv[3]);
  double c_const = convertTo<double>(argv[4]);
  double g_const = convertTo<double>(argv[5]);
  double t_const = convertTo<double>(argv[6]);

  IntensityReader           input_intensity_file(input_intensity_filename);
  IntensityFilter_Constant  constfilter(a_const,c_const,g_const,t_const);
  input_intensity_file.reopen();

  IntensityWriter output_intensity_file(output_intensity_filename);

  for(;!input_intensity_file.eof();) {

    // An intensitysequence is a single line of the intensity filename
    IntensitySequence is = input_intensity_file.get_next();
    
    /// *** FITLERS GO HERE ***
    is = constfilter.filter(is);

    if(!input_intensity_file.eof()) {
        output_intensity_file.write(is);
    }
  }
 
  return 0;
}
