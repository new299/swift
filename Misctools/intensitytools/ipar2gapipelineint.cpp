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

#include "intensity.h"
#include "iparintensityreader.h"
#include "intensitywriter.h"
#include "stringify.h"

using namespace std;

int main(int argc,char **argv) {

  if(argc < 3) {
    cerr << "iparintensity2gapipelineint <ipar intensity file> <ipar position file> <output gapipeline intfile>" << endl;
    return 0; 
  }

  const char* input_ipar_intensity_filename  = argv[1];
  const char* input_ipar_positions_filename  = argv[2];
  const char* output_intfile_filename        = argv[3];

  string iparfile = argv[1];
  size_t position = iparfile.find_last_of("/");

  string inputbasename = iparfile.substr(position+1,iparfile.size()-position-1);

  cout << "inputbasename: " << inputbasename << endl;
  
  int lane=0;
  int tile=0;

  try {
    lane = convertTo<int>(inputbasename.substr(2,1));
    tile = convertTo<int>(inputbasename.substr(4,4));
  } catch (std::exception e) {
  }

  IparIntensityReader input_intensity_file(input_ipar_intensity_filename,input_ipar_positions_filename,lane,tile);

  IntensityWriter out(output_intfile_filename);
  for(;!input_intensity_file.eof();) {
    // An intensitysequence is a single line of the intensity filename
    IntensitySequence is = input_intensity_file.get_next();

    if(!input_intensity_file.eof()) {
      out.write(is);
    }
  }
 
  return 0;
}
