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

#ifndef INTENSITYTOOLS_INTENSITY2PAIRWISEPLOT_H
#define INTENSITYTOOLS_INTENSITY2PAIRWISEPLOT_H

#include "intensityreader.h"
#include <vector>
#include <iostream>

using namespace std;

void intensity2pairwiseplot(char *intensityfile,char *plotprefix) {
  // Intensity file reader
  IntensityReader intensity_file(intensityfile);
 
  // open plot files
  ofstream plot_file_cg((string(plotprefix) + "cg").c_str());
  ofstream plot_file_ca((string(plotprefix) + "ca").c_str());
  ofstream plot_file_ct((string(plotprefix) + "ct").c_str());
  ofstream plot_file_ga((string(plotprefix) + "ga").c_str());
  ofstream plot_file_gt((string(plotprefix) + "gt").c_str());
  ofstream plot_file_at((string(plotprefix) + "at").c_str());

  for(;!intensity_file.eof();) {
   
    // An intensitysequence is a single line of the intensity file
    IntensitySequence is = intensity_file.get_next();
    
    if(!intensity_file.eof()) {
      for(IntensitySequence::intensities_type::iterator i = is.intensities.begin();i != is.intensities.end();i++) {

        // dump intensities to plots in order to create a scatter graph (with gnuplot)
        plot_file_cg << (*i).getbase(ReadIntensity::base_c) << " " << (*i).getbase(ReadIntensity::base_g) << endl;
        plot_file_ca << (*i).getbase(ReadIntensity::base_c) << " " << (*i).getbase(ReadIntensity::base_a) << endl;
        plot_file_ct << (*i).getbase(ReadIntensity::base_c) << " " << (*i).getbase(ReadIntensity::base_t) << endl;
        plot_file_ga << (*i).getbase(ReadIntensity::base_g) << " " << (*i).getbase(ReadIntensity::base_a) << endl;
        plot_file_gt << (*i).getbase(ReadIntensity::base_g) << " " << (*i).getbase(ReadIntensity::base_t) << endl;
        plot_file_at << (*i).getbase(ReadIntensity::base_a) << " " << (*i).getbase(ReadIntensity::base_t) << endl;
      }
    }
  }
 
  // close plots
  plot_file_cg.close();
  plot_file_ca.close();
  plot_file_ct.close();
  plot_file_ga.close();
  plot_file_gt.close();
  plot_file_at.close();
}

#endif
