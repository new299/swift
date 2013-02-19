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

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "FastaReader.h"
#include "FastaSequence.h"
#include "SmallAlign.h"
#include "stringify.h"
#include "CommandLine.h"
#include "./FastqReader.h"
#include "./FastqWriter.h"
#include "./FastPrbReader.h"
#include "./ProbabilitySequence.h"

using namespace std;
  
int main(int argc,char **argv) {
 
  CommandLine *parms = CommandLine::Instance();


  parms->add_valid_parm("mapfile","Mapping between qualities");
  parms->add_valid_parm("in","Input file");
  parms->add_valid_parm("out","Output file");

  cout << parms->usage() << endl;

  parms->process_args(argc,argv);

  // Load map
  cout << "opening map: " << parms->get_parm("mapfile") << endl;
  ifstream mapfile(parms->get_parm("mapfile").c_str());

  if(!mapfile.is_open()) return 1;
  map<int,int> qualitymap;
  for(;!mapfile.eof();) {
    
    string qual1;
    mapfile >> qual1;

    string qual2;
    mapfile >> qual2;

    if(qual1.compare("") != 0) {
      cout << "qual1: " << qual1 << "    qual2: " << qual2 << endl;
      qualitymap[convertTo<int>(qual1)] = convertTo<int>(qual2);
    }
  }

  FastqReader<BaseProbability<4> > fastq_in(parms->get_parm("in"));
  FastqWriter<BaseProbability<4> > fastq_out(parms->get_parm("out"));

  if((!fastq_in.is_open()) || (!(fastq_out.is_open()))) {
    cerr << "Unable to open fastq files" << endl;
    return 1;
  }

  for(;!fastq_in.eof();) {
    ProbabilitySequence<BaseProbability<4> > p = fastq_in.get_sequence();
    string qualitystr = fastq_in.get_last_quality_string();

    string new_qualitystr;
    for(int n=0;n<qualitystr.length();n++) {
      int q = qualitystr[n];
      q -= fastq_in.get_quality_offset();

      int newq = qualitymap[q];

      newq += fastq_out.get_quality_offset();
      
      new_qualitystr.push_back(newq);
    }

    fastq_out.write(p.get_id(),p.string_sequence(),new_qualitystr);
  }

  fastq_in.close();
  fastq_out.close();
}
