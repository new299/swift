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

#ifndef SWIFT_SMALLALIGN_H
#define SWIFT_SMALLALIGN_H

#include "SequenceAlignment.h"
#include "sequence_utils.h"

template<class sequence_type=string,class _prec=int>
class SmallAlign {

public:
  SmallAlign(sequence_type reference_in,
             bool add_reverse,
             bool circular_genome_in)
           : circular_genome(circular_genome_in) {

    total_ref_size = 0;

    if(reference_in.length() < 6000) {
      references.push_back(reference_in);
      if(add_reverse) {
        references.push_back(sequence_reverse(reference_in));
      }

    } 

    if(all_reference_size() > 40000) {
      cerr << "Error: SmallAlign is design to process SMALL sequences (6000bp and less) reference size exceeds this, alignment will not be performed." << endl;
    }
  }
  
  SmallAlign(bool circular_genome_in)
           : circular_genome(circular_genome_in) {
  }

  /// no bounds checking
  _prec get_offset(_prec referencenum) {
    _prec offset=0;
    for(int n=0;n < referencenum;n++) {
      offset += references[n].size();
    }

    return offset;
  }

  bool add_reference(sequence_type reference,bool add_reverse) {

    if(reference.length() > 6000) {
      cerr << "Error: SmallAlign is design to process SMALL sequences (6000bp and less) this reference will not be used." << endl;
      references.clear();
      return true;
    }

    references.push_back(reference);

    if(add_reverse) {
      references.push_back(sequence_reverse(reference));
    }

    if(all_reference_size() > 40000) {
      cerr << "Error: Total reference size exceeds 40000bp, SmallAlign is not designed to process sequences of this size, alignment will not be performed" << endl;
      references.clear();
    }

    return true;
  }

  _prec get_num_references() const {
    return references.size();
  }
  
  _prec all_reference_size() const {

    _prec total=0;
    for(typename vector<sequence_type>::const_iterator i=references.begin();i != references.end();i++) {
      total += (*i).size();
    }

    return total;
  }

  string get_ref(_prec reference,_prec position,_prec length) {
    string s;
    
    if(references[reference].size() == 0) return s; // Zeroed size reference will break this code

    for(int n=0,p=0;n<length;n++,p++) {
      if(static_cast<unsigned int>(position+p) >= references[reference].size()) {position=0; p=0;}
      s.push_back(references[reference][position+p]);
    }

    return s;
  }

  SequenceAlignment<_prec> align(sequence_type sequence) {

    _prec best_score           = -1;
    _prec best_score_reference = -1;
    _prec best_score_position  = -1;
    bool  best_unique          = false;
    vector<bool> best_score_matchstring;

    for(typename vector<sequence_type>::const_iterator i = references.begin();i != references.end();i++) {
      for(typename sequence_type::const_iterator j = (*i).begin();j != (*i).end();j++) {
        // Compare sequence to this position
        vector<bool> matchstring;
        _prec score = comparison(sequence.begin(),sequence.end(),(*i).begin(),j,(*i).end(),matchstring,(*i).size()-best_score);

        if(score >= best_score) {
          if(score==best_score) best_unique=false; else best_unique=true;
          best_score             = score;
          best_score_reference   = i-references.begin();
          best_score_position    = j-(*i).begin();
          best_score_matchstring = matchstring;
        }
      }
    }

    return SequenceAlignment<_prec>(best_score_reference,best_score_position,best_score,best_score_matchstring,best_unique);
  }

  _prec comparison(typename sequence_type::const_iterator i_begin,
                   typename sequence_type::const_iterator i_end,
                   typename sequence_type::const_iterator j_start, ///< Start of the reference sequence of which j is part, used for circular genomes
                   typename sequence_type::const_iterator j_begin,
                   typename sequence_type::const_iterator j_end,
                   vector<bool> &matchstring,
                   int max_errs) {
    
    typename sequence_type::const_iterator i_cur = i_begin;
    typename sequence_type::const_iterator j_cur = j_begin;
    
    _prec score=0;
    int errs=0;
    matchstring.clear();
    for(;(i_cur != i_end) && (circular_genome || (j_cur != j_end));i_cur++) {
      if(circular_genome) if(j_cur == j_end) j_cur=j_start;
      if(j_cur != j_end) {          // Copes with 0 length sequence
        if((*i_cur) == (*j_cur)) {
          score++;
          matchstring.push_back(true);
          if(errs > max_errs) return -1; 
        } else {errs++; matchstring.push_back(false);}
      }
      if(j_cur != j_end) j_cur++;          // Iterate, it not at the end already... Copes with 0 length sequence
    }
    
    // no penalty for running off the end of the sequence, in non-circular genome

    return score;
  }

private:
  vector<sequence_type> references;
  size_t total_ref_size;
  bool circular_genome;
};

#endif
