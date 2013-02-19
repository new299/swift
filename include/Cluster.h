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

#ifndef SWIFT_CLUSTER_H
#define SWIFT_CLUSTER_H

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include "stringify.h"
#include "ClusterPosition.h"
#include "ReadIntensity.h"
#include "ProbabilitySequence.h"
#include <algorithm>
  
// To use alternate allocator:
/*
  typedef std::vector<intensities_type,__gnu_cxx::new_allocator<intensities_type> >        signal_vec_type;  ///< Type which holds a vector of signal intensities.
  typedef std::vector<intensities_type,__gnu_cxx::new_allocator<intensities_type> >        noise_vec_type;   ///< Type which holds a vector of noise intensities.
  map<string,signal_vec_type,less<string>,__gnu_cxx::new_allocator<pair<string,signal_vec_type>  > > signal_vec;                 ///< Map of Vectors of Signal values
  map<string,noise_vec_type ,less<string>,__gnu_cxx::new_allocator<pair<string,signal_vec_type>  > >  noise_vec;                  ///< Map of Vectors of noise estimates
*/

template <class _prec = float, class _position_prec = int> 
class Cluster {

public:
  typedef ReadIntensity<_prec>                 intensities_type; ///< Type which holds 4 intensities (one for each base).
  typedef std::vector<intensities_type>        signal_vec_type;  ///< Type which holds a vector of signal intensities.
  typedef std::vector<intensities_type>        noise_vec_type;   ///< Type which holds a vector of noise intensities.
  typedef ProbabilitySequence<>                sequence_type;    ///< Type which holds bases called

  bool valid;                                             ///< Valid or not, not sure if I'm happy with this being public.

private:
  string last_process_signal;                             ///< This is the string identifer for the last signal   processing operation that occured.
  string last_process_sequence;                           ///< This is the string identifer for the last sequence processing operation that occured.
  ClusterPosition<_position_prec> m_position;             ///< The position of this cluster
  map<string,signal_vec_type,less<string> > signal_vec;   ///< Map of Vectors of Signal values
  map<string,noise_vec_type ,less<string> > noise_vec;    ///< Map of Vectors of noise estimates
  vector<string>              signal_ids;                 ///< List of signal ids
  vector<string>              sequence_ids;               ///< List of signal ids
  map<string,sequence_type>   sequences;                  ///< Called sequences

public:
  Cluster() : valid(true),
              last_process_signal("RAW"),
              last_process_sequence("BASECALL") { // last_process should probably left blank here and set when the sequence/signals are inserted.
    signal_vec["RAW"] = signal_vec_type();
    noise_vec["RAW"]  = noise_vec_type();
    signal_ids.push_back("RAW");
    // sequence_ids.push_back("BASECALL"); // Should I be adding this here, or let BaseCaller do it?
  }
  
  inline bool is_valid() const {
    return valid;
  }

  inline void set_valid(bool v) {
    valid = v;
  }

  bool add_signal(string identifier) {
    
    if(signal_vec.find(identifier) == signal_vec.end()) {
      signal_vec[identifier] = signal_vec_type();
      noise_vec [identifier] = noise_vec_type();
      last_process_signal = identifier;
      signal_ids.push_back(identifier);
    }

    return true;
  }

  bool add_sequence(string identifier) {
    
    if(sequences.find(identifier) == sequences.end()) {
      sequences[identifier] = sequence_type();
      last_process_sequence = identifier;
      sequence_ids.push_back(identifier);
    }

    return true;
  }

  inline void delete_signal(const string &identifier) {
   
    signal_vec_type mst;
    noise_vec_type  mnt;

    if(signal_vec.find(identifier) == signal_vec.end()) cerr << "Error in Cluster.h: identifier not found during delete: " << identifier << endl;
    signal(identifier).swap(mst);
    noise (identifier).swap(mnt);

    signal_vec.erase(identifier);
    noise_vec .erase(identifier);
    
    if(signal_vec.find(identifier) != signal_vec.end()) cerr << "Error in Cluster.h: identifier not deleted in delete: " << identifier << endl;
    
    std::remove(signal_ids.begin(),signal_ids.end(),identifier);
  }

  /// Accessor for signal, no bounds checking
  signal_vec_type &signal(const string &identifier) {
    return (*signal_vec.find(identifier)).second;
  }

  /// Accessor for noise, no bounds checking
  noise_vec_type &noise(const string &identifier) {
    return (*noise_vec.find(identifier)).second;
  }
  
  /// Const Accessor for signal, no bounds checking
  const signal_vec_type &const_signal(const string &identifier) const {

    if(signal_vec.find(identifier) == signal_vec.end()) {
      cerr << "ERROR TAG: " << identifier << " NOT FOUND IN CLUSTER" << endl;
    }

    return (*signal_vec.find(identifier)).second;
  }

  /// Const Accessor for noise, no bounds checking
  const noise_vec_type &const_noise(const string &identifier) const {
    return (*noise_vec.find(identifier)).second;
  }
  
  /// Accessor for sequence, no bounds checking
  sequence_type &sequence(string identifier) {
    return (*sequences.find(identifier)).second;
  }
  
  /// Accessor for sequence, no bounds checking
  sequence_type &sequence(string identifier, string tag) {
    sequence_type &s = (*sequences.find(identifier)).second;
    s.set_id(tag + string(":") + stringify(m_position.x) + string(":") + stringify(m_position.y));
  
    return s;
  }

  /// Const Accessor for signal, no bounds checking
  const sequence_type &const_sequence(string identifier,string tag="") const {
    return (*sequences.find(identifier)).second;
  }

  /// Accessor for RAW signal, no bounds checking
  signal_vec_type &raw_signal() {
    return (*signal_vec.find("RAW")).second;
  }

  /// Accessor for RAW noise, no bounds checking
  noise_vec_type &raw_noise() {
    return (*noise_vec.find("RAW")).second;
  }
  
  /// Accessor for RAW signal, no bounds checking, const version
  const signal_vec_type &const_raw_signal() const {
    return (*signal_vec("RAW")).second;
  }

  /// Accessor for RAW noise, no bounds checking, const version
  const noise_vec_type &const_raw_noise() const {
    return (*noise_vec("RAW")).second;
  }

  // Accessor for LAST PROCESSED signal
  signal_vec_type &processed_signal() {
    return (*signal_vec(last_process_signal)).second;
  }
  
  // Accessor for LAST_PROCESSED noise
  noise_vec_type &processed_noise() {
    return (*signal_vec(last_process_signal)).second;
  }
  
  // Const accessor for LAST_PROCESSED signal
  const signal_vec_type &const_processed_signal() const {
    return (*signal_vec(last_process_signal)).second;
  }

  // Const accessor for LAST_PROCESSED noise
  const noise_vec_type &const_processed_noise() const {
    return (*noise_vec(last_process_signal)).second;
  }

  const string last_processed_signal() const {
    return last_process_signal;
  }
  
  const string last_processed_sequence() const {
    return last_process_sequence;
  }

  const vector<string> get_signal_ids() const {
    return signal_ids;
  }
  
  const vector<string> get_sequence_ids() const {
    return sequence_ids;
  }

  const ClusterPosition<_position_prec> get_position() const {
    return m_position;
  }

  void set_position(const ClusterPosition<_position_prec> &newpos) {
    m_position = newpos;
  }
 

  /// Compare similarity of two cluster signals
  /// If they are within some similarity threshold, return true
  /// First attempt is similarity based on base with maximum
  /// intensity.
  int similarity(const Cluster<_prec> &other,string signalid) {

    const signal_vec_type &mysignal    =       const_signal(signalid);
    const signal_vec_type &othersignal = other.const_signal(signalid);

    if(mysignal.size() != othersignal.size()) return false;

    int similar_bases=0;
    for(size_t n=0;n<mysignal.size();n++) {
      if(n < othersignal.size()) {
        if(mysignal[n].max_base() == othersignal[n].max_base()) {
          similar_bases++;
        }
      }
    }
    
    return similar_bases;
  }


  // Return minimum Purity of bases between start_base and end_base.
  // No bounds checking.
  _prec min_purity(int start_base,int end_base,string signalid) const {
    
    if(end_base < start_base) return 0;
    _prec min_purity = const_signal(signalid)[start_base].purity();
    for(int i=start_base+1;i <= end_base;i++) {
      if(const_signal(signalid)[i].purity() < min_purity) min_purity = const_signal(signalid)[i].purity();
    }
    return min_purity;
  }
  
  // Return second lowest purity
  _prec min_purity2(int start_base,int end_base,string signalid) const {
    
    _prec min_purity  = const_signal(signalid)[start_base].purity();
    _prec min_purity2 = min_purity;
    for(int i=start_base+1;i <= end_base;i++) {

      _prec current_purity = const_signal(signalid)[i].purity();

      if(current_purity <= min_purity) {
        min_purity2 = min_purity;
        min_purity = current_purity;
      } else

      if(current_purity <= min_purity2) {
        min_purity2 = current_purity;
      }
    }
    return min_purity2;
  }

  bool min_purity_greaterthaneq(int start_base,int end_base,string signalid,_prec threshold) const {
    for(int i=start_base;i <= end_base;i++) {
      if(const_signal(signalid)[i].purity() < threshold) return false;
    }
    return true;
  }

  _prec average_purity(string signalid="RAW",int length=-1) const {

    _prec purity_sum=0;

    size_t n=0;
    for(typename signal_vec_type::const_iterator i = const_signal(signalid).begin();(i != const_signal(signalid).end()) && (n < length);i++,n++) {
      purity_sum += (*i).purity();
    }

    return purity_sum/const_signal(signalid).size();
  }

  bool off_edge(int threshold,string signalid) const {
    int offedgecount=0;
    for(typename signal_vec_type::const_iterator i = const_signal(signalid).begin();i != const_signal(signalid).end();i++) {
      if((*i).off_edge()) offedgecount++;
    }
  
    if(offedgecount > threshold) return true; else return false;
  }
  
  bool any_off_edge(string signalid, int threshold) const {
    int offedgecount=0;
    for(typename signal_vec_type::const_iterator i = const_signal(signalid).begin();i != const_signal(signalid).end();i++) {
      if((*i).any_off_edge()) offedgecount++;
    }
  
    if(offedgecount > threshold) return true; else return false;
  }
  
  bool first_off_edge(string signalid,unsigned int howmany) const {
    int offedgecount=0;

    //TODO: Use an iterator here
    for(size_t i = 0;i < howmany;i++) {
      if(const_signal(signalid)[i].off_edge()) offedgecount++;
    }
  
    if(offedgecount > 0) return true;
                 else    return false;
  }

  _prec average_peaksignal(string signalid) const {
    _prec sum=0;
    for(typename signal_vec_type::const_iterator i = const_signal(signalid).begin();i != const_signal(signalid).end();i++) {
      sum += (*i).max_intensity();
    }

    return sum/signal_vec.size();
  }
  
  _prec average_peaksignal_within(string signalid,int bases=0) const {
    _prec sum=0;
    for(int n=0;n < bases;n++) {
      sum += const_signal(signalid)[n].max_intensity();
    }

    return sum/bases;
  }


  string dump_gapipelinestr(string signalid) const {
    string s;
    s += "1 1 ";
    s += get_position().as_string();
    s += " ";
    for(typename signal_vec_type::const_iterator i = const_signal(signalid).begin();i != const_signal(signalid).end();i++) {
      s += (*i).as_string();
      s += " ";
    }

    return s;
  }
  
  bool read_gapipelinestr(string signalid,string line) {
    stringstream in(line);
    
    string s;
    
    // discard lane and tile number
    in >> s;
    in >> s;

    in >> s;
    int x = convertTo<int>(s);
    
    in >> s;
    int y = convertTo<int>(s);
    set_position(ClusterPosition<_position_prec>(x,y));
    
    signal(signalid).clear();
    for(;!in.eof();) {
      string sa;
      in >> sa;

      string sc;
      in >> sc;

      string sg;
      in >> sg;

      bool end_of_line = in.eof();
      string st;
      in >> st;
     
      if(!end_of_line) {
        double a = convertTo<double>(sa);
        double c = convertTo<double>(sc);
        double g = convertTo<double>(sg);
        double t = convertTo<double>(st);

        bool a_offedge=false;
        bool c_offedge=false;
        bool g_offedge=false;
        bool t_offedge=false;

        if(a == 0.0) a_offedge=true;
        if(c == 0.0) c_offedge=true;
        if(g == 0.0) g_offedge=true;
        if(t == 0.0) t_offedge=true;
        
        signal(signalid).push_back(ReadIntensity<_prec>(a,c,g,t,a_offedge,c_offedge,g_offedge,t_offedge));
      }
    }

    return true;
  }

};
 
#include "Cluster.cpp"

#endif
