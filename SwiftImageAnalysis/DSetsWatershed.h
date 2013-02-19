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

#ifndef SWIFTIMAGEANALYSIS_DSETSWATERSHED
#define SWIFTIMAGEANALYSIS_DSETSWATERSHED

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

// This is a modified disjoint set algorithm which allows multiple parents
// and therefore multiple canonical ids to be associated with a set.

class DSetsWatershed {
public:
  vector<vector<int> > parents;

   DSetsWatershed(int dsets_size,bool compact_in=true) : do_compact(compact_in) {
     parents.clear();
     parents.insert(parents.begin(),dsets_size,vector<int>());
   }

   void makeparent(int c,int l) {
     // cout << "makeparent called: " << c << "," << l << endl;
     parents[c].push_back(l);
   }

   void add(int p) {
     
     // cout << "parents.size(): " << parents[p].size() << "[0]=";
     if(parents[p].size() == 0) {
        parents[p].push_back(p);
     } //else {
       //cout << parents[p][0] << endl;
    // }
   }

   vector<int> getparents(int c) {
     vector<int> r;
     //cout << "get parents called: " << c << endl;
     for(vector<int>::iterator i = parents[c].begin();i != parents[c].end();i++) {
       if((*i) != c) {
         vector<int> p = getparents((*i));
         if(p.size() == 0) r.push_back(*i); // if you can resolve the parent, leave it as is.
                      else r.insert(r.begin(),p.begin(),p.end());
       } else {
         r.clear();
         r.push_back(c);
         return r;
       }
     }

     parents[c] = r;
     if(do_compact) compact(c);

     return parents[c];
   }

   // Removes duplicates from parents list.
   // Is this required? does it cripple the algorithm?
   void compact(int c) {
     //cout << "compact called: " << c << endl; 
     vector<int> pnew;
     sort(parents[c].begin(),parents[c].end());

     if(parents[c].size() != 0)
     for(unsigned int n=0;n<(parents[c].size()-1);n++) {
       if(parents[c][n] != parents[c][n+1]) {
         pnew.push_back(parents[c][n]);
       }
     }

     if(parents[c].size() != 0) pnew.push_back(parents[c][parents[c].size()-1]);
   
     parents[c]=pnew;
   }

   void dump(ostream &out) {
     for(unsigned int n=0;n<parents.size();n++) {
       out << n << ": ";
       for(unsigned int na=0;na<parents[n].size();na++) {
         out << parents[n][na] << ",";
       }
       out << endl;
     }
   }

private:
   bool do_compact;
};

#endif
