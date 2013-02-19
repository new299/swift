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

#include "intensityfilter_flowcell.h"
#include "stringify.h"

using namespace std;

int main(int argc,char **argv) {

  if(argc < 12) {
    cout << "intensity2flowcellremove <input intensity file> <input noise file> <input positions file> <output good intensity file> <output good noise file> <output good positions file> <output bad intensity file> <output bad noise file> <output bad positions file> <num bins> <threshold>" << endl;
    cout << "Performs deepmagic to remove flow cell walls." << endl;
    cout << "ok, what it actually does is bin the image in the x direction (similar to intensity2xybins)" << endl;
    cout << "it then detects a hump in these bins and removes all these intensities from the intensity file." << endl;
    cout << "The ``hump'' is determined by the average, if the bin contain more than average*threshold clusters is a high bin" << endl;
    cout << "if it contains less, it's a low bin. To detect a hump the bins must go from low to high to low again (starts in state low)" << endl;
    cout << "within the first 25% of the image." << endl;
    return 0; 
  }

  const char* input_intensity_filename  = argv[1];
  const char* input_noise_filename      = argv[2];
  const char* input_positions_filename   = argv[3];
  const char* output_good_intensity_filename = argv[4];
  const char* output_good_noise_filename     = argv[5];
  const char* output_good_positions_filename  = argv[6];
  const char* output_bad_intensity_filename = argv[7];
  const char* output_bad_noise_filename     = argv[8];
  const char* output_bad_positions_filename  = argv[9];

  int    flowcellfilter_numbins            = convertTo<int>(argv[10]);
  double flowcellfilter_threshold          = convertTo<double>(argv[11]);

  cout << "Input intensity file: " << input_intensity_filename << endl;
  cout << "Input noise file    : " << input_noise_filename     << endl;
  cout << "Input position file : " << input_positions_filename  << endl;
  

  IntensityReader input_intensity_file(input_intensity_filename);
  IntensityFilter_Flowcell flowcellfilter(input_intensity_file,flowcellfilter_numbins,flowcellfilter_threshold,cout);

  input_intensity_file.reopen();
  ifstream input_noise_file    (input_noise_filename);
  ifstream input_positions_file(input_positions_filename);

  IntensityWriter output_good_intensity_file(output_good_intensity_filename);
  ofstream        output_good_noise_file    (output_good_noise_filename);
  ofstream        output_good_positions_file(output_good_positions_filename);

  IntensityWriter output_bad_intensity_file;
  ofstream        output_bad_noise_file;
  ofstream        output_bad_positions_file;

  int discarded=0;
  for(;!input_intensity_file.eof();) {

    // An intensitysequence is a single line of the intensity filename
    IntensitySequence is = input_intensity_file.get_next();
    
    string noise_line;
    getline(input_noise_file,noise_line);

    string positions_line;
    getline(input_positions_file,positions_line);

    /// *** FITLERS GO HERE ***
    is = flowcellfilter.filter(is);

    if(!input_intensity_file.eof()) {
      if(is.valid) {
        output_good_intensity_file.write(is);
        output_good_noise_file     << noise_line << endl;
        output_good_positions_file << positions_line << endl;
      } else {
        
        if(discarded==0) {
          output_bad_intensity_file.open(output_bad_intensity_filename);
          output_bad_noise_file    .open(output_bad_noise_filename);
          output_bad_positions_file.open(output_bad_positions_filename);
        }

        output_bad_intensity_file.write(is);
        output_bad_noise_file      << noise_line << endl;
        output_bad_positions_file << positions_line << endl;
        
        discarded++;
      }
    }
  }
 
  cout << "Discarded: " << discarded << endl;
  return 0;
}
