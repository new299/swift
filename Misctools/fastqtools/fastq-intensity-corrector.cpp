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

#include "fastqreader.h"
#include "fastqwriter.h"
#include "../intensitytools/intensityreader.h"
#include "../intensitytools/intensitywriter.h"

#include <string>

using namespace std;

int main(int argc,char **argv) {

  if(argc < 2) {
    cout << "fastqcount <fastq file> <intensity input file> <intensity output file>" << endl;
   return 1;
  }

  FastqReader fastq_in(argv[1]); 
  const char *intensity_input_filename = argv[2];
  const char *intensity_output_filename = argv[3];


  fastq_in.open();

  cout << "Input file  : " << argv[1] << endl;

  bool eof=false;
 

  vector<int> a_count(500,0); // hardcoded bad!
  vector<int> t_count(500,0); // hardcoded bad!
  vector<int> g_count(500,0); // hardcoded bad!
  vector<int> c_count(500,0); // hardcoded bad!
  vector<int> n_count(500,0); // hardcoded bad!
  int max_size=0;
  for(;!eof;) {
    ScoredSequence s = fastq_in.next_line(eof);
    
    for(int n=0;(n < s.sequence.size());n++) {
      if(s.sequence[n] == 'A') a_count[n]++; else
      if(s.sequence[n] == 'T') t_count[n]++; else
      if(s.sequence[n] == 'G') g_count[n]++; else
      if(s.sequence[n] == 'C') c_count[n]++; else
      n_count[n]++;
    }

    if(s.sequence.size() > max_size) max_size = s.sequence.size();
  }
  fastq_in.close();
 
  cout << "RAW" << endl;
  cout << "A count:"; 
  for(int i = 0;i < max_size;i++) cout << a_count[i]  << " ";
  cout << endl;
  
  cout << "T count:";
  for(int i = 0;i < max_size;i++) cout << t_count[i]  << " ";
  cout << endl;
  
  cout << "G count:";
  for(int i = 0;i < max_size;i++) cout << g_count[i]  << " ";
  cout << endl;
  
  cout << "C count:";
  for(int i = 0;i < max_size;i++) cout << c_count[i]  << " ";
  cout << endl;
  
  cout << "N count:";
  for(int i = 0;i < max_size;i++) cout << n_count[i]  << " ";
  cout << endl;
  cout << endl;

  vector<int> total(500,0);
  for(int i=0;i<max_size;i++) {
    total[i] = a_count[i] + t_count[i] + g_count[i] + c_count[i] + n_count[i];
  }

  vector<double> a_count_normalised(500,0);
  vector<double> t_count_normalised(500,0);
  vector<double> g_count_normalised(500,0);
  vector<double> c_count_normalised(500,0);
  for(int i=0;i<max_size;i++) {
    a_count_normalised[i] = static_cast<double>(a_count[i])/static_cast<double>(total[i]);
    t_count_normalised[i] = static_cast<double>(t_count[i])/static_cast<double>(total[i]);
    g_count_normalised[i] = static_cast<double>(g_count[i])/static_cast<double>(total[i]);
    c_count_normalised[i] = static_cast<double>(c_count[i])/static_cast<double>(total[i]);
  }


  cout << "Normalised" << endl;
  cout << "A count:"; 
  for(int i = 0;i < max_size;i++) cout << a_count_normalised[i]  << " ";
  cout << endl;
  
  cout << "T count:"; 
  for(int i = 0;i < max_size;i++) cout << t_count_normalised[i]  << " ";
  cout << endl;

  cout << "G count:"; 
  for(int i = 0;i < max_size;i++) cout << g_count_normalised[i]  << " ";
  cout << endl;
  
  cout << "C count:"; 
  for(int i = 0;i < max_size;i++) cout << c_count_normalised[i]  << " ";
  cout << endl;
  
  vector<double> at_correct(500,0);
  vector<double> gc_correct(500,0);
  
  for(int i=0;i<max_size;i++) {
    at_correct[i] = (((double)a_count[i]/(double)total[i])+((double)t_count[i]/(double)total[i]))/2;
    gc_correct[i] = (((double)g_count[i]/(double)total[i])+((double)c_count[i]/(double)total[i]))/2;
  }

  cout << "CORRECT TO" << endl;
  cout << "A/T correct: "; 
  for(int i = 0;i < max_size;i++) cout << at_correct[i] << " ";
  cout << endl;
  
  cout << "G/C correct: ";
  for(int i = 0;i < max_size;i++) cout << gc_correct[i] << " ";
  cout << endl;
 
  vector<double> a_correction(500,0);
  vector<double> t_correction(500,0);
  vector<double> g_correction(500,0);
  vector<double> c_correction(500,0);
  
  for(int i=0;i<max_size;i++) {
    a_correction[i] = at_correct[i]/a_count_normalised[i];
    t_correction[i] = at_correct[i]/t_count_normalised[i];
    g_correction[i] = gc_correct[i]/g_count_normalised[i];
    c_correction[i] = gc_correct[i]/c_count_normalised[i];
  }

  cout << "Correction: " << endl;
  cout << "A correction:";
  for(int i = 0;i < max_size;i++) cout << a_correction[i]  << " ";
  cout << endl;
  
  cout << "T correction:";
  for(int i = 0;i < max_size;i++) cout << t_correction[i]  << " ";
  cout << endl;
  
  cout << "G correction:";
  for(int i = 0;i < max_size;i++) cout << g_correction[i]  << " ";
  cout << endl;
  
  cout << "C correction:";
  for(int i = 0;i < max_size;i++) cout << c_correction[i]  << " ";
  cout << endl;


  // Intensity file reader
  IntensityReader input_intensity_file(intensity_input_filename);
  IntensityWriter output_intensity_file(intensity_output_filename);

  for(;!input_intensity_file.eof();) {

    // An intensitysequence is a single line of the intensity file
    IntensitySequence is = input_intensity_file.get_next();

    if(!input_intensity_file.eof()) {
      for(int i=0;i< is.intensities.size();i++) {
        is.intensities[i].bases[ReadIntensity::base_a] = is.intensities[i].bases[ReadIntensity::base_a]*a_correction[i]; 
        is.intensities[i].bases[ReadIntensity::base_t] = is.intensities[i].bases[ReadIntensity::base_t]*t_correction[i]; 
        is.intensities[i].bases[ReadIntensity::base_g] = is.intensities[i].bases[ReadIntensity::base_g]*g_correction[i]; 
        is.intensities[i].bases[ReadIntensity::base_c] = is.intensities[i].bases[ReadIntensity::base_c]*c_correction[i]; 
      }

      output_intensity_file.write(is);
    }
  }

  input_intensity_file.close();
  output_intensity_file.close();
}
