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
#include "./FastPrbReader.h"
#include "./ProbabilitySequence.h"

using namespace std;
  
int main(int argc,char **argv) {
 
  CommandLine *parms = CommandLine::Instance();

  parms->add_valid_parm("ref","Reference sequence for alignment");
  parms->add_valid_parm("seq","Sequences to align (fastq)");
  parms->add_valid_parm("table","Table mapping between qualities");
  parms->add_valid_parm("snplist","A list of SNP positions to exclude during table generation");

  if(argc < 2) {
    cout << parms->usage() << endl; 
    return 1;
  }

  parms->process_args (argc, argv);

  SmallAlign<> aligner(true);

  //Perform alignment, calculate error rate
  
  //Load SNP table

  vector<int> snp_positions;
  ifstream snpfile(parms->get_parm("snplist").c_str());

  if(snpfile.is_open()) {
    string snppos;
    snpfile >> snppos;

    int snpval = convertTo<int>(snppos);

    snp_positions.push_back(snpval);
  }


  // Load reference sequence
  FastaReader<> reference_fastafile(parms->get_parm("ref"));
  bool openok = reference_fastafile.open();

  if(openok == false) cerr << "Could not open reference file: " << parms->get_parm("ref") << endl;
  else                cerr << "Opened reference file: " << parms->get_parm("ref") << endl;
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
  
  FastqReader<BaseProbability<4> > sequences_fastq(parms->get_parm("seq"));

  openok=true;
  if(openok == false) cerr << "Could not open sequences file: " << parms->get_parm("seq") << endl;
  else                cerr << "Opened sequences file: " << parms->get_parm("seq") << endl;
  fasta_eof=false;
 

  int max_cycles=40;
  int max_quality=100;
  int printevery=1000;
  vector<int> errors_by_cycle(max_cycles,0);
  vector<int> errors_by_quality(max_quality,0);
  vector<int> quality_count(max_quality,0);
  vector<int> total_by_cycle (max_cycles,0);
  
  vector<int> errors_count(max_cycles,0);

  if(openok) {
    int n=0;
    for(;!sequences_fastq.eof();) {
      ProbabilitySequence<BaseProbability<4> > seq = sequences_fastq.get_sequence();
      
      SequenceAlignment<> alm;
      alm = aligner.align(seq.string_sequence());

      if(alm.score > 30) 
      for(unsigned int m=0;m<alm.matchstring.size();m++) {
        
        // slow...
        bool this_snp = false;
        for(int spos=0;spos<snp_positions.size();spos++) {
          if((alm.contig == 0) && (alm.position+m) == snp_positions[spos]) {
            this_snp = true;
          }
          if((alm.contig == 1) && (alm.position+m) == (aligner.get_references()[1].size()-snp_positions[spos])) {
            this_snp = true;
          }
        }

        if(this_snp == true) cout << "snp pos!" << endl;
        if(this_snp == false) {
          if(alm.matchstring[m] == false) {
            errors_by_quality[sequences_fastq.get_last_quality_string()[m]-sequences_fastq.get_quality_offset()]++;
            errors_by_cycle[m]++;
          }
          quality_count[sequences_fastq.get_last_quality_string()[m]-sequences_fastq.get_quality_offset()]++;
          total_by_cycle[m]++;
        }
      }
      errors_count[alm.score]++;

      if((n%printevery == 0) || sequences_fastq.eof()) {
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
        
        // Errors by quality
        cout << "Error by quality" << endl;
        for(int quality=0;(quality < max_quality);quality++) {
          if(quality_count[quality] != 0) {
            double p_err = (static_cast<double>(errors_by_quality[quality])/static_cast<double>(quality_count[quality]));
            cout << setw(5) << quality << " ";
            cout << setw(5) << errors_by_quality[quality] << " ";
            cout << setw(5) << quality_count[quality] << " ";
            cout << setw(10) << p_err << " ";
            cout << setw(5) << -10*log10(p_err) << endl;
          }  
        }
      }
      n++;
    }
  } else {
    cerr << "Could not open sequences file" << endl;
  }

  ofstream ftable(parms->get_parm("table").c_str());

  for(int quality=0;(quality < max_quality);quality++) {
    if(quality_count[quality] != 0) {
      double p_err = (static_cast<double>(errors_by_quality[quality])/static_cast<double>(quality_count[quality]));
      ftable << quality;
      ftable  << " ";
      ftable <<  static_cast<int>((-10*log10(p_err))+0.5) << endl;
    }
  }
 
  ftable.close();

}
