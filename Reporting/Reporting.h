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

#ifndef SWIFT_REPORTING_H
#define SWIFT_REPORTING_H

#include <iostream>
#include <fstream>
#include <vector>
#include "Cluster.h"
#include "clusterfilter_purity.h"
#include "clusterfilter_purityavg.h"
#include "SmallAlign.h"
#include "FastaReader.h"
#include "Memstats.h"
#include "Timetagger.h"

using namespace std;

template<class _prec=double>
class Reporting {
public:

  Reporting(vector<Cluster<_prec> > &clusters,    ///< Generate report on these clusters
            bool   do_alignment=false,            ///< True if you want an a alignment
            string reference_genome_filename="",  ///< Reference genome for this tile
            bool   is_paired=false,               ///< Is this a paired end run?
            int    pair_break_pos=0,              ///< The pair break position
            size_t align_every=100,               ///< Align every Nth read
            string tiletag_in="",                 ///< Tag to write at top of output
            string contams_filename="",           ///< Contaminants reference for this tile
            ostream &err_in=std::cerr             ///< Write debug/errors here
           ) : m_clusters(clusters),
               m_reference_genome_filename(reference_genome_filename),
               m_is_paired(is_paired),
               m_pair_break_pos(pair_break_pos),
               m_align_every(align_every),
               m_tiletag(tiletag_in),
               m_contams_filename(contams_filename),
               m_do_alignment(do_alignment),
               err(err_in) {
    generate_stats();
  }

  bool generate_stats() {
    
    // Crosstalk plots 1st cycle, last cycle

    // Average intensities, all and called
    Memstats mem ("calc_average_intensities");
    calc_average_intensities("FINAL");

    // Error rates, PF/Non-PF/Cycle, actually also calculate number of pf/non-pf bases
    mem.start ("calc_error_rates");
    calc_error_rates();
    mem.stop();
  
    return true;
  }

  bool calc_average_intensities(string signal_id) {

    average_all_intensities.clear();
    average_all_intensities.insert(average_all_intensities.begin(),m_clusters[0].const_signal(signal_id).size(),vector<_prec>(ReadIntensity<_prec>::base_count,0));
    
    average_called_intensities.clear();
    average_called_intensities.insert(average_called_intensities.begin(),m_clusters[0].const_signal(signal_id).size(),vector<_prec>(ReadIntensity<_prec>::base_count,0));

    //TODO: really don't like calculating number of clusters like this
    for(size_t cycle=0;cycle < m_clusters[0].const_signal("FINAL").size();cycle++) {
      
      vector<_prec> total_all_signal(ReadIntensity<_prec>::base_count,0);
      vector<int>   total_all_count(ReadIntensity<_prec>::base_count,0);
      vector<_prec> total_called_signal(ReadIntensity<_prec>::base_count,0);
      vector<int> total_called_count(ReadIntensity<_prec>::base_count,0);
      
      for(typename vector<Cluster<_prec> >::const_iterator i=m_clusters.begin();i != m_clusters.end();i++) {
        for(int base=0;base < ReadIntensity<_prec>::base_count;base++) {
          total_all_signal[base] += (*i).const_signal(signal_id)[cycle].get_base(base);
          total_all_count[base]++;

          typename ReadIntensity<_prec>::base_type maxb = (*i).const_signal(signal_id)[cycle].max_base();
          _prec maxb_val = (*i).const_signal(signal_id)[cycle].get_base((*i).const_signal(signal_id)[cycle].max_base());

          total_called_signal[maxb] += maxb_val;
          total_called_count[maxb]++;
        }
      }

      for(int base=0;base < ReadIntensity<_prec>::base_count;base++) {
        average_all_intensities   [cycle][base] = total_all_signal[base]/total_all_count[base];
        average_called_intensities[cycle][base] = total_called_signal[base]/total_called_count[base];
      }
     }
  
    return true;
  }

  ///TODO: Currently not using contams file
  bool calc_error_rates() {
    SmallAlign<> aligner(true);

    pf_total_bases=0;
    nonpf_total_bases=0;
    
    pf_total_bases_aligned=0;
    nonpf_total_bases_aligned=0;

    pf_total_reads=0;
    nonpf_total_reads=0;
    
    pf_total_reads_aligned=0;
    nonpf_total_reads_aligned=0;

    pf_total_score=0;
    nonpf_total_score=0;

    pf_unique=0;
    nonpf_unique=0;

    pf_score_by_cycle.clear();
    nonpf_score_by_cycle.clear();

    //TODO: again getting cycle number from first cluster and assuming they are all the same is bad
    pf_score_by_cycle.insert(pf_score_by_cycle.begin(),m_clusters[0].const_sequence("FINAL").size(),0);
    nonpf_score_by_cycle.insert(nonpf_score_by_cycle.begin(),m_clusters[0].const_sequence("FINAL").size(),0);

    if(m_do_alignment) {
      // Load reference sequence
      FastaReader<> fastafile(m_reference_genome_filename);
      bool openok = fastafile.open();

      if(openok == false) err << "Could not open reference file: " << m_reference_genome_filename << endl;
      else                err << "Opened reference file: " << m_reference_genome_filename << endl;
      bool fasta_eof=false;

      if(openok)
      for(;!fasta_eof;) {
        FastaSequence reference = fastafile.next_sequence(fasta_eof);
        if(!fasta_eof) {
          aligner.add_reference(reference.sequence,true);
        }

        if(reference.sequence.size() > 6000) {
          m_do_alignment=false;
          err << "SmallAlign, Swift's built in aligner is a brute force aligner, and should not be used to align sequencers greater than 6Kb (at most)" << endl;
          err << "No alignment will be performed" << endl;
        }
      }

      if(aligner.all_reference_size() <= 0) m_do_alignment = false;
    }

    Memstats mem ("Alignments");
    size_t align_count = 0;
    for(typename vector<Cluster<_prec> >::const_iterator i = m_clusters.begin();i != m_clusters.end();i++) {

      string seq;
       
      int ends=0;
      if(m_is_paired) ends=2; else ends=1;

      for(int end=0;end<ends;end++) {

        int offset=0;
        if(m_is_paired) {
          if(end == 0) {
            seq = (*i).const_sequence("FINAL").trim(0,m_pair_break_pos-1).get_sequence_string();
            offset=0;
          }

          if(end == 1) {
            seq = (*i).const_sequence("FINAL").trim(m_pair_break_pos,-1 ).get_sequence_string();
            offset = m_pair_break_pos;
          }
        } else {
          seq = (*i).const_sequence("FINAL").get_sequence_string();
          offset=0;
        }
      
        if((*i).is_valid()) {
          pf_total_bases += seq.length();
          pf_total_reads++;
        } else {
          nonpf_total_bases += seq.length();
          nonpf_total_reads++;
        }
      
        if((align_count % m_align_every) == 0) {

          if(m_do_alignment) {
            SequenceAlignment<> alm;
            alm = aligner.align(seq);
            if((*i).is_valid()) {
              pf_total_score += alm.score;
            
              for(unsigned int m=0;m<alm.matchstring.size();m++) {
                if(alm.matchstring[m] == false) pf_score_by_cycle[m+offset]++;
              }
              if(alm.unique == true) pf_unique++;
            

              pf_total_bases_aligned += seq.length();
              pf_total_reads_aligned++;

            } else {
              
              if(alm.unique == true) nonpf_unique++;
              nonpf_total_score += alm.score;
              
              nonpf_total_bases_aligned += seq.length();
              nonpf_total_reads_aligned++;
              
              for(unsigned int m=0;m<alm.matchstring.size();m++) {
                if(alm.matchstring[m] == false) nonpf_score_by_cycle[m+offset]++;
              }
            }
          }
        }

      }
      align_count++;
    }
    mem.stop();
  
    return true;
    
  }

  bool write_report_human(ostream &out) {
    out << "PF Reads   : " << pf_total_reads << " (both ends)" << endl;
    out << "NonPF Reads: " << nonpf_total_reads << " (both ends)" << endl;
    
    if(m_do_alignment) {
      out << "Error rates based on aligning every " << m_align_every << "th read" << endl;
      out << "PF Error rate: " << ((static_cast<double>(pf_total_bases_aligned)-static_cast<double>(pf_total_score))/static_cast<double>(pf_total_bases_aligned))*100 << "%" << endl;

      out << "PF Errors by Cycle" << endl;
      for(unsigned int cycle=0;cycle < pf_score_by_cycle.size();cycle++) {
        out << cycle << " " << ((static_cast<double>(pf_score_by_cycle[cycle]))/(static_cast<double>(pf_total_reads_aligned)))*100 << endl;
      }
    }
    out << endl;
  
    return true;
  }

  bool write_report_file(string filename,string extraxml) {
  
    ofstream out(filename.c_str());

    out << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl;
    out << "<runreport>" << endl;
    out << "<dataset>SwiftReport</dataset>" << endl;
    out << "<report tag=\"" << m_tiletag << "\" version=\"$Id: Reporting.h 182 2009-04-29 14:37:12Z nw3 $\" />" << endl;

    out << "<summary>" << endl;

    out << "<pfbases     value=\"" << pf_total_bases    << "\" />" << endl;
    out << "<nonpfbases  value=\"" << nonpf_total_bases << "\" />" << endl;

    out << "<pfreads     value=\"" << pf_total_reads    << "\" />" << endl;
    out << "<nonpfreads  value=\"" << nonpf_total_reads << "\" />" << endl;
    

    out << "</summary>" << endl;
    if(m_do_alignment) {
      out << "<alignment>" << endl;
    
      // PF Data
      out << "<pf>" << endl;
      out << "<goodbases      value=\"" << pf_total_score    << "\" />" << endl;
      out << "<uniquereads    value=\"" << pf_unique         << "\" />" << endl;

      out << "<errorbycycle>" << endl;
      for(unsigned int cycle=0;cycle < pf_score_by_cycle.size();cycle++) {
        out << "<row cycle=\"" << cycle+1 << "\" errors=\"" << pf_score_by_cycle[cycle] << "\" />" << endl;
      }
      out << "</errorbycycle>" << endl;
      out << "</pf>" << endl;

      // NonPF Data
      out << "<nonpf>" << endl;
      out << "<goodbases      value=\"" << nonpf_total_score    << "\" />" << endl;
      out << "<uniquereads    value=\"" << nonpf_unique         << "\" />" << endl;

      out << "<errorbycycle>" << endl;
      for(unsigned int cycle=0;cycle < nonpf_score_by_cycle.size();cycle++) {
        out << "<row cycle=\"" << cycle+1 << "\" errors=\"" << nonpf_score_by_cycle[cycle] << "\" />" << endl;
      }
      out << "</errorbycycle>" << endl;
      out << "</nonpf>" << endl;

      out << "</alignment>" << endl;
    }

    out << "<intensities>" << endl;
    
    out << "<averageall>" << endl;
    for(unsigned int cycle=0;cycle < average_all_intensities.size();cycle++) {
      out << "<row cycle=\"" << cycle+1 << "\"";
      for(int base=0;base<ReadIntensity<_prec>::base_count;base++) {
        out << " " << ReadIntensity<_prec>::base_name[base] << "=\"" << average_all_intensities[cycle][base] << "\"";
      }
      out << " />" << endl;
    }
    out << "</averageall>" << endl;
    
    out << "<averagecalled>" << endl;
    for(unsigned int cycle=0;cycle < average_called_intensities.size();cycle++) {
      out << "<row cycle=\"" << cycle+1 << "\"";
      for(int base=0;base<ReadIntensity<_prec>::base_count;base++) {
        out << " " << ReadIntensity<_prec>::base_name[base] << "=\"" << average_called_intensities[cycle][base] << "\"";
      }
      out << " />" << endl;
    }
    out << "</averagecalled>" << endl;

    out << "</intensities>" << endl;
    out << extraxml << endl;        // any other xml that's passed to this method
    out << "</runreport>" << endl;
    return true;
  }

private:
  const vector<Cluster<_prec> > &m_clusters;
  string      m_reference_genome_filename;
  bool        m_is_paired;
  int         m_pair_break_pos;
  size_t      m_align_every;
  string      m_tiletag;
  string      m_contams_filename;


  bool m_do_alignment;

  vector<_prec> ctalk_ac;
  vector<_prec> ctalk_ag;
  vector<_prec> ctalk_at;
  vector<_prec> ctalk_cg;
  vector<_prec> ctalk_ct;
  vector<_prec> ctalk_gt;

  vector<vector<_prec> > average_all_intensities;
  vector<vector<_prec> > average_called_intensities;

  vector<int> nonpf_score_by_cycle;
  int         nonpf_total_bases;
  int         nonpf_total_bases_aligned;
  int         nonpf_total_score;
  int         nonpf_total_reads;
  int         nonpf_total_reads_aligned;
  int         nonpf_unique;

  vector<int> pf_score_by_cycle;
  int         pf_total_bases;
  int         pf_total_bases_aligned;
  int         pf_total_score;
  int         pf_total_reads;
  int         pf_total_reads_aligned;
  int         pf_unique;


  ostream &err;                               ///< Error output will be writen here, set to cerr in constructor default
};

// This needs to be included here because this is a template class
// #include "Reporting.cpp"

#endif
