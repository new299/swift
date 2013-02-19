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
#include "FastqWriter.h"
#include "FastaSequence.h"
#include "SmallAlign.h"
#include "stringify.h"
#include "CommandLine.h"

using namespace std;
  
int main(int argc,char **argv) {
  
  CommandLine parms;

  parms.add_valid_parm("ref","Reference sequence for alignment");
  parms.add_valid_parm("seq","Sequences to align (fasta)");

  cout << "SmallAlign" << endl;
  cout << parms.usage() << endl; 

  parms.process_args(argc, argv);

  SmallAlign<> aligner(true);

  //Perform alignment, calculate error rate
  
  // Load reference sequence
  FastaReader<> reference_fastafile(parms.get_parm("ref"));
  bool openok = reference_fastafile.open();

  if(openok == false) cerr << "Could not open reference file: " << parms.get_parm("ref") << endl;
  else                cerr << "Opened reference file: " << parms.get_parm("ref") << endl;
  bool fasta_eof=false;

  if(openok) {
    for(;!fasta_eof;) {
      FastaSequence reference = reference_fastafile.next_sequence(fasta_eof);
      if(!fasta_eof) {
        // Forward and Reverse? (true/false)
        aligner.add_reference(reference.sequence,true);
      }
    }
  }
  
  FastaReader<> sequences_fasta(parms.get_parm("seq"));
  openok = sequences_fasta.open();

  if(openok == false) cerr << "Could not open sequences file: " << parms.get_parm("seq") << endl;
  else                cerr << "Opened sequences file: " << parms.get_parm("seq") << endl;
  fasta_eof=false;

  int max_cycles=40;
  int printevery=1000;
  vector<int> errors_by_cycle(max_cycles,0);
  vector<int> total_by_cycle (max_cycles,0);
  
  vector<int> errors_count(max_cycles,0);


  size_t align_limit = 1000;

  if(openok) {
    int n=0;
    for(;!fasta_eof;) {
      FastaSequence seq = sequences_fasta.next_sequence(fasta_eof);
      if(!fasta_eof) {
        SequenceAlignment<> alm;
        alm = aligner.align(seq.sequence);

        for(unsigned int m=0;m<alm.matchstring.size();m++) {
          if(alm.matchstring[m] == false) errors_by_cycle[m]++;
          total_by_cycle[m]++;
        }
        errors_count[alm.score]++;
      }

      if((n%printevery == 0) || fasta_eof) {
        // Errors by cycle
        cout << "Error rates" << endl;
        for(int cycle=0;(cycle < max_cycles) && (total_by_cycle[cycle] > 0);cycle++) {
          cout << setw(5) << cycle << " ";
          cout << setw(10) << (static_cast<double>(errors_by_cycle[cycle])/static_cast<double>(total_by_cycle[cycle]))*100 << " ";
          cout << setw(10) << errors_by_cycle[cycle] << " ";
          cout << setw(10) << total_by_cycle[cycle] << endl;
        }
        cout << "Score count" << endl;
        for(int errs=0;(errs < max_cycles);errs++) {
          cout << setw(5) << errs << " ";
          cout << setw(10) << errors_count[errs] << endl;
        }
      }
      n++;

      if(n > align_limit) fasta_eof = true;
    }
  } else {
    cerr << "Could not open sequences file" << endl;
  }

}
