////////////////////////////////////////////////////////////////////////////////////////
//
// A C++ interface to gnuplot.
//
//
// The interface uses pipes and so wont run on a system that doesn't have POSIX pipe support
// Tested on Windows and Linux
//
// Version history:
// 0. C interface written by N. Devillard (27/01/03)
//    http://ndevilla.free.fr/gnuplot/
// 1. C++ interface: direct translation from the C interface by Rajarshi Guha (07/03/03)
//    http://cheminfo.informatics.indiana.edu/~rguha/code/cc++/
// 2. corrections for Win32 compatibility by V.Chyzhdzenka (20/05/03)
//    chyzhdzenka@mail.ru
// 3. corrections for Win32 and Linux compatibility, some member functions added
//    by M. Burgis (22/01/08)
//
// Requirements:
// * gnuplot has to be installed (http://www.gnuplot.info/download.html)
// * set Path-Variable for Gnuplot path (e.g. C:/program files/gnuplot/bin; /usr/bin)
//     or set Gnuplot path with: Gnuplot::set_GNUPlotPath(const std::string &path)
//
////////////////////////////////////////////////////////////////////////////////////////


#ifndef _GNUPLOT_PIPES_H_
#define _GNUPLOT_PIPES_H_


#include <string>
#include <vector>
#include <stdexcept>  // for std::runtime_error class in GnuplotException
#include <cstdio>     // for FILE


//declare classes in global namespace


class GnuplotException : public std::runtime_error
{
    public:
        GnuplotException(const std::string &msg) : std::runtime_error(msg){}
};



class Gnuplot
{
private:

  // member data
  FILE                    *gnucmd;
  bool                     valid;        // validation of gnuplot session
  bool                     two_dim;      // true = 2d, false = 3d
  int                      nplots;       // number of plots in session
  std::string              pstyle;       // functions and data are displayed in a defined styles
  std::vector<std::string> tmpfile_list; // list of created tmpfiles

  // static data
  static int               tmpfile_num;        // number of all tmpfiles (number of tmpfiles restricted)
  static std::string       m_sGNUPlotFileName; // Name of executed GNUPlot file
  static std::string       m_sGNUPlotPath;     // gnuplot path

  // member functions (auxiliary functions)
  void           init(bool checkForX11);              // get_program_path(); and popen();
  std::string    create_tmpfile(std::ofstream &tmp);  // creates tmpfile and returns its name

  //static functions
  static bool    get_program_path();                                   // gnuplot path found?
  static bool    file_exists(const std::string &filename, int mode=0); // checks if file exists


public:

  // optional: set Gnuplot path manual
  //   for windows: path with slash "/" not backslash "\"
  static bool set_GNUPlotPath(const std::string &path);

  Gnuplot(bool checkForX11);
  
  // set a style during construction
  Gnuplot(const std::string &style = "points", 
	  bool checkForX11 = true);

  // plot a single std::vector at one go
  Gnuplot(const std::vector<double> &x,
	  const std::string &title = "",
	  const std::string &style = "points",
	  const std::string &labelx = "x",
	  const std::string &labely = "y", 
	  bool checkForX11 = true);

  // plot pairs std::vector at one go
  Gnuplot(const std::vector<double> &x,
	  const std::vector<double> &y,
	  const std::string &title = "",
	  const std::string &style = "points",
	  const std::string &labelx = "x",
	  const std::string &labely = "y", 
	  bool checkForX11 = true);

  // plot triples std::vector at one go
  Gnuplot(const std::vector<double> &x,
	  const std::vector<double> &y,
	  const std::vector<double> &z,
	  const std::string &title = "",
	  const std::string &style = "points",
	  const std::string &labelx = "x",
	  const std::string &labely = "y",
	  const std::string &labelz = "z", 
	  bool checkForX11 = true);

  // destructor: needed to delete temporary files
  ~Gnuplot();

  // send a command to gnuplot
  Gnuplot& cmd(const std::string &cmdstr);
  Gnuplot& operator<<(const std::string &cmdstr);

  // sets terminal type to windows / x11
  Gnuplot& showonscreen();

  // saves a gnuplot session to a postscript file, filename without extension
  Gnuplot& savetops(const std::string &filename = "gnuplot_output");

  // set line style (some of these styles require additional information)
  // lines, points, linespoints, impulses, dots, steps, fsteps, histeps,
  // errorbars, labels, xerrorbars, yerrorbars, xyerrorbars, errorlines,
  // xerrorlines, yerrorlines, xyerrorlines, boxes, histograms, filledcurves,
  // boxerrorbars, boxxyerrorbars, financebars, candlesticks, vectors, image,
  // rgbimage, pm3d.
  Gnuplot& set_style(const std::string &stylestr = "points");

  // scales the size of the points used in plots
  Gnuplot& set_pointsize(const double pointsize = 1.0);

  // turns grid on/off
  Gnuplot& set_grid();
  Gnuplot& unset_grid();

  // turns on/off log scaling for the specified axes
  Gnuplot& set_xlogscale(const double base = 10);
  Gnuplot& set_ylogscale(const double base = 10);
  Gnuplot& set_zlogscale(const double base = 10);
  Gnuplot& unset_xlogscale();
  Gnuplot& unset_ylogscale();
  Gnuplot& unset_zlogscale();

  // set sampling rate of functions, or for interpolating data
  Gnuplot& set_samples(const int samples = 100);
  // set isoline density (grid) for plotting functions as surfaces (for 3d plots)
  Gnuplot& set_isosamples(const int isolines = 10);

  // enables/disables hidden line removal for surface plotting (for 3d plot)
  Gnuplot& set_hidden3d();
  Gnuplot& unset_hidden3d();

  // switches legend on/off
  // position: inside/outside, left/center/right, top/center/bottom, nobox/box
  Gnuplot& set_legend(const std::string &position = "default");
  Gnuplot& unset_legend();

  // sets and clears the title of a gnuplot session
  Gnuplot& set_title(const std::string &title = "");
  Gnuplot& unset_title();

  // set y and x axis labels
  Gnuplot& set_ylabel(const std::string &label = "x");
  Gnuplot& set_xlabel(const std::string &label = "y");
  Gnuplot& set_zlabel(const std::string &label = "z");

  // set axis - ranges
  Gnuplot& set_xrange(const int iFrom,
		      const int iTo);
  Gnuplot& set_yrange(const int iFrom,
		      const int iTo);
  Gnuplot& set_zrange(const int iFrom,
		      const int iTo);

  // set palette range
  Gnuplot& set_cbrange(const int iFrom, const int iTo);


  // plot a single std::vector: x
  //   from file
  Gnuplot& plotfile_x(const std::string &filename,
		      const std::string &title = "");
  //   from std::vector
  Gnuplot& plot_x(const std::vector<double> &x,
		  const std::string &title = "");


  // plot x,y pairs: x y
  //   from file
  Gnuplot& plotfile_xy(const std::string &filename,
		       const std::string &title = "");
  //   from std::vector
  Gnuplot& plot_xy(const std::vector<double> &x,
		   const std::vector<double> &y,
		   const std::string &title = "");


  // plot x,y pairs with dy errorbars: x y dy
  //   from file
  Gnuplot& plotfile_xy_err(const std::string &filename,
			   const std::string &title = "");
  //   from std::vector
  Gnuplot& plot_xy_err(const std::vector<double> &x,
		       const std::vector<double> &y,
		       const std::vector<double> &dy,
		       const std::string &title = "");


  // plot x,y,z triples: x y z
  //   from file
  Gnuplot& plotfile_xyz(const std::string &filename,
			const std::string &title = "");
  //   from std::vector
  Gnuplot& plot_xyz(const std::vector<double> &x,
		    const std::vector<double> &y,
		    const std::vector<double> &z,
		    const std::string &title = "");


  // plot an equation of the form: y = ax + b, you supply a and b
  Gnuplot& plot_slope(const double a,
		      const double b,
		      const std::string &title = "");

  // plot an equation supplied as a std::string y=f(x), write only the function f(x) not y=
  // functions: abs(x), abs(x), acos(x), acosh(x), arg(x), asin(x), asinh(x), atan(x),
  //  atan2(y,x), atanh(x), besj0(x), besj1(x), besy0(x), besy1(x), ceil(x), cos(x), cosh(x),
  //  erf(x), erfc(x), exp(x), floor(x), gamma(x), ibeta(p,q,x), inverf(x), igamma(a,x),
  //  imag(x), invnorm(x), int(x), lambertw(x), lgamma(x), log(x), log10(x), norm(x), rand(x),
  //  real(x), sgn(x), sin(x), sinh(x), sqrt(x), tan(x), tanh(x),
  Gnuplot& plot_equation(const std::string &equation,
			 const std::string &title = "");

  // plot an equation supplied as a std::string z=f(x,y), write only the function f(x,z) not z=
  Gnuplot& plot_equation3d(const std::string &equation,
			   const std::string &title = "");


  // plot image
  Gnuplot& plot_image(const unsigned char * ucPicBuf,
		      const int iWidth,
		      const int iHeight,
		      const std::string &title = "");


  // resets a gnuplot session (next plot will erase previous ones)
  Gnuplot& reset_plot();

  // validation of gnuplot session
  bool is_valid();

};


#endif
