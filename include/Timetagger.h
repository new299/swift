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

#ifndef TIMETAGGER_H_
#define TIMETAGGER_H_

#include <sys/time.h>
#include "stringify.h"

// Return a timestamp string containing the current time to millisecond granularity.
// Typically, caller would create a Timetagger object once, then repeatedly invoke
// its methods to fetch timestamps.

using namespace std;

class Timetagger {
public:

  Timetagger () {
    str();            // put something meaningful in member vars
  }

  int now () {
    gettimeofday(&m_tim, NULL);
    return m_tim.tv_sec;
  }

  string str () {

    string m_timestamp;

    tm tm_time;
    
    now();
    localtime_r(&m_tim.tv_sec, &tm_time);
    
    m_timestamp += stringify(tm_time.tm_hour,3);
    m_timestamp += ":";
    m_timestamp += stringify(tm_time.tm_min,2);
    m_timestamp += ":";
    m_timestamp += stringify(tm_time.tm_sec,3);
    m_timestamp += ":";
    m_timestamp += stringify(m_tim.tv_usec/1000,3);   
    m_timestamp += " ";

    // sprintf (m_timestamp, "%02d:%02d:%02d.%03ld ",
    //    tm_time.tm_hour, tm_time.tm_min, tm_time.tm_sec, m_tim.tv_usec/1000);
    
    return m_timestamp;
    
  }

private:

  timeval m_tim;
  char    m_timestamp[32];

};

#endif /*TIMETAGGER_H_*/
