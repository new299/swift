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

#ifndef MEMSTATS_H_
#define MEMSTATS_H_

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include "stringify.h"
#include "Timetagger.h"

using namespace std;

// Interrogate /proc/<pid>/statm to extract the amount of memory user by this process.
// TODO: statm contains 7 values, describing usage of various sorts of memory. Make 
// them all available.

class Memstats {
public:

  // Create object, don't start the clock, don't print anything.
  Memstats () : m_start_time(0), m_start_mem(0), m_message("") {
    int pid = getpid();
    m_proc_name = "/proc/" + stringify(pid) + "/statm";
  }
  
  // Create object, start the clock, print start message.
  Memstats (string message) : m_start_time(0), m_message(message) {
    int pid = getpid();
    m_proc_name = "/proc/" + stringify(pid) + "/statm";
    start(message);
  }
  
  ~Memstats () {        // print stop message when object goes out of scope
    stop();
  }
  
  void start (string message = "") {

    if (m_start_time > 0) {      // if clock is running, stop it and print message
      stop();
    }
    
    m_message = message; 
    m_start_time = m_tt.now();
    m_start_mem = get_total();
    
    cout << m_tt.str() << "Memstats start:       " 
      << right << setw(8) << m_start_mem << "            "    // need to match spacing of stop message
      << m_message << endl;
  
  }
  
  void stop () {
    
    if (m_start_time > 0) {
      
      int stop_mem = get_total();
      int del_mem  = stop_mem - m_start_mem;
      int del_time = m_tt.now() - m_start_time;
      
      cout << m_tt.str() << "Memstats stop: (" << right << setw(4) << del_time << ") "
        << right << setw(8) << stop_mem 
        << " (" << right << setw(8) << del_mem << ") "  
        << m_message << endl;
      
    }
    
    m_start_time = 0;              // clock is not running
    
  }
  
	int get_total () {
		ifstream proc(m_proc_name.c_str());
		int total = 0;
		proc >> total;     // memory in pages
		total *= 4;        // convert to Kbytes  TODO: dynamically determine pagesize
		return total;
	}
	
	void print (string message = "") {
		cout << m_tt.str() << "Memstats (Kbytes): " << right << setw(8) << get_total() << " " << message << endl;
	}

private:
  
  int        m_start_time;
  int        m_start_mem;
	string     m_proc_name;
	string     m_message;
  Timetagger m_tt;
	
};

#endif /*MEMSTATS_H_*/
