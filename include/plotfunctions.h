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

#ifndef SWIFT_PLOTFUNCTIONS
#define SWIFT_PLOTFUNCTIONS

#include <iostream>
#include <fstream>
#include <vector>
#include "CrossTalkCorrection.h"
#include "gnuplot_i.hpp"
#include "Cluster.h"
#include "plotfunctions.h"

using namespace std;
  
void crosstalk_plot(const vector<Cluster<> > &clusters,          ///< Clusters to plot
                    string signalid,                             ///< Which signal id to plot
                    int cycle,ReadIntensity<>::base_type base_x, ///< X Axis base for crosstalk plot
                    ReadIntensity<>::base_type base_y,           ///< Y Axis base for crosstalk plot
                    string plottitle) {                          ///< Plot title/data name  
  Gnuplot *g = new Gnuplot(plottitle);
  Gnuplot &crosstalkplot = *g;
  vector<double> plot_x;
  vector<double> plot_y;

  for(vector<Cluster<> >::const_iterator i=clusters.begin();i != clusters.end();i++) {
    if((*i).valid) {
      plot_x.push_back((*i).const_signal(signalid)[cycle].get_base(base_x));
      plot_y.push_back((*i).const_signal(signalid)[cycle].get_base(base_y));
    }
  }
  
  crosstalkplot.set_grid();
  crosstalkplot.reset_plot();
  crosstalkplot.set_style("points");
  crosstalkplot.plot_xy(plot_x,plot_y,plottitle + ReadIntensity<>::base_name[base_x] + "/" + ReadIntensity<>::base_name[base_y]);
  
  // Gnuplot's distructor deletes the plot temp file before Gnuplot has a chance to plot the graph
  // so for the moment I'm not deleting Gnuplot, this sucks and should be fixed.
}

void plotxy(const vector<double> &plot_x,
            const vector<double> &plot_y,
            string plottitle) {                          ///< Plot title/data name  
  Gnuplot *g = new Gnuplot(plottitle);
  Gnuplot &crosstalkplot = *g;

  crosstalkplot.set_grid();
  crosstalkplot.reset_plot();
  crosstalkplot.set_style("points");
  crosstalkplot.plot_xy(plot_x,plot_y,plottitle);
  
  // Gnuplot's distructor deletes the plot temp file before Gnuplot has a chance to plot the graph
  // so for the moment I'm not deleting Gnuplot, this sucks and should be fixed.

}

#endif
