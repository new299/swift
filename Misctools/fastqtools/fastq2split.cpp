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
#include "stringify.h"
#include <map>

#include <string>

using namespace std;

string extract_tag(string str,string delim,int tag_num) {
  
  size_t pos=0;

  for(size_t n=0;n<tag_num;n++) {
    pos = str.find(delim,pos+1);
  }
  size_t e_pos = str.find(delim,pos+1);

  if(pos == string::npos) return "";
  string o = str.substr(pos+1,e_pos-pos-1);
 
  return o;
}

int main(int argc,char **argv) {

  if(argc < 3) {
    cout << "fast2split <fastq file> <delim> <tag number> <output base name>" << endl;
   return 0;
  }

  ofstream outfile(argv[2]);

  cout << "Input file  : " << argv[1] << endl;
  cout << "Output file : " << argv[2] << endl;

  string delim = argv[2];
  int tag_num = convertTo<int>(argv[3]);
  string outbase = argv[4];

  // 1. Find all tag values
  FastqReader fastq_in (argv[1]); 
  fastq_in.open();

  map<string,int> tag_map;
  bool eof=false;
  for(;!eof;) {
    ScoredSequence s = fastq_in.next_line(eof);
   
    string tag = extract_tag(s.id,delim,tag_num);
    tag_map[tag]++;
  }

  ofstream tagcount((outbase + string(".count")).c_str());
  for(map<string,int>::const_iterator i=tag_map.begin();i != tag_map.end();i++) {
    tagcount << (*i).first << " " << (*i).second << endl;
  }
  tagcount.close();

  // 2. Create tag output files
  map<string,FastqWriter *> output_files;
  for(map<string,int>::const_iterator i=tag_map.begin();i != tag_map.end();i++) {
    output_files[(*i).first] = new FastqWriter(outbase + string(".") + (*i).first);
    output_files[(*i).first]->open();
  }

  // 3. parse input file for tags
  FastqReader fastq_in2 (argv[1]); 
  fastq_in2.open();
  eof = false;
  for(;!eof;) {
    ScoredSequence s = fastq_in2.next_line(eof);
   
    if(!eof) {
      string tag = extract_tag(s.id,delim,tag_num);
      output_files[tag]->write(s);
    }
  }

  // 4. close output files
  for(map<string,int>::const_iterator i=tag_map.begin();i != tag_map.end();i++) {
    output_files[(*i).first]->close();
  }

  fastq_in.close();
  fastq_in2.close();
}













