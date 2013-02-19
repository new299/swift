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

#include "intensity2xybins.h"
#include "stringify.h"

using namespace std;

int main(int argc,char **argv) {

  if(argc < 2) {
    cout << "intensity2xybins <input intensity file> <num bins>" << endl;
    cout << "Count number of reads in strips across tile in x and y direction <num bins> is the number of strips" << endl;
    return 0; 
  }

  vector<int> xbins;
  vector<int> ybins;

  double binsize=0;
  IntensityReader input_intensity_file(argv[1]);
  intensity2xybins(input_intensity_file,convertTo<int>(argv[2]),binsize,xbins,ybins);

  cout << "Bin size: " << binsize << endl;
  cout << "XBINS: " << endl;
  for(vector<int>::iterator i = xbins.begin();i != xbins.end();i++) cout << (*i) << " ";
  cout << endl;
  
  cout << "YBINS: " << endl;
  for(vector<int>::iterator i = ybins.begin();i != ybins.end();i++) cout << (*i) << " ";
  cout << endl;

  return 0;
}
