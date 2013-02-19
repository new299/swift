////////////////////////////////////////////
//
// A C++ interface to gnuplot.
//
// This is a direct translation from the C interface
// written by N. Devillard (which is available from
// http://ndevilla.free.fr/gnuplot/).
//
// As in the C interface this uses pipes and so wont
// run on a system that doesn't have POSIX pipe
// support
//
// Rajarshi Guha
// <rajarshi@presidency.com>
//
// 07/03/03
//
////////////////////////////////////////////
//
// A little correction for Win32 compatibility
// and MS VC 6.0 done by V.Chyzhdzenka
//
// Notes:
// 1. Added private method Gnuplot::init().
// 2. Temporary file is created in th current
//    folder but not in /tmp.
// 3. Added #indef WIN32 e.t.c. where is needed.
// 4. Added private member m_sGNUPlotFileName is
//    a name of executed GNUPlot file.
//
// Viktor Chyzhdzenka
// e-mail: chyzhdzenka@mail.ru
//
// 20/05/03
//
////////////////////////////////////////////
//
// corrections for Win32 and Linux compatibility
//
// some member functions added:
//  set_GNUPlotPath
//  create_tmpfile, get_program_path, file_exists
//  operator<<, plotfile_*, plot_xy_err
//  plot_equation3d, savetops, showonscreen
// set, unset: pointsize, grid, *logscale, samples,
//  isosamples, hidden3d, legend, title, cbrange
//
// Markus Burgis
//
// 22/01/08
//
////////////////////////////////////////////


#include <fstream>              // for std::ifstream
#include <sstream>              // for std::ostringstream
#include <list>                 // for std::list
#include <cstdio>               // for FILE, fputs(), fflush(), popen()
#include <cstdlib>              // for getenv()
#include "gnuplot_i.hpp"

/*
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__) //defined for 32 and 64-bit environments
 #include <io.h>                // for _access(), _mktemp()
 #define GP_MAX_TMP_FILES    27 // 27 temporary files it's Microsoft restriction
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)//all UNIX-like OSs (Linux, *BSD, MacOSX, Solaris, ...)
 #include <unistd.h>            // for access(), mkstemp()
 #define GP_MAX_TMP_FILES    64
#else
 #error unsupported or unknown operating system
#endif
*/

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__) //defined for 32 and 64-bit environments
 #include <io.h>                // for _access(), _mktemp()
 #define GP_MAX_TMP_FILES    27 // 27 temporary files it's Microsoft restriction
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__) //all UNIX-like OSs (Linux, *BSD, MacOSX, Solaris, ...)
 #include <unistd.h>            // for access(), mkstemp()
 #define GP_MAX_TMP_FILES    64
#else
 #error unsupported or unknown operating system
#endif



//----------------------------------------------------------------------------------
//
// initialize static data
//
int Gnuplot::tmpfile_num = 0;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
std::string Gnuplot::m_sGNUPlotFileName = "pgnuplot.exe";
std::string Gnuplot::m_sGNUPlotPath = "C:/program files/gnuplot/bin/";
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
std::string Gnuplot::m_sGNUPlotFileName = "gnuplot";
std::string Gnuplot::m_sGNUPlotPath = "/usr/bin/";
#endif


//----------------------------------------------------------------------------------
//
// define static member function
//
bool Gnuplot::set_GNUPlotPath(const std::string &path)
{

  std::string tmp = path + "/" + Gnuplot::m_sGNUPlotFileName;


#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
  if ( Gnuplot::file_exists(tmp,0) ) // check existence
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    if ( Gnuplot::file_exists(tmp,1) ) // check existence and execution permission
#endif
      {
        Gnuplot::m_sGNUPlotPath = path;
        return true;
      }
    else
      {
        Gnuplot::m_sGNUPlotPath.clear();
        return false;
      }
}


//----------------------------------------------------------------------------------
//
// A string tokenizer taken from http://www.sunsite.ualberta.ca/Documentation/
// /Gnu/libstdc++-2.90.8/html/21_strings/stringtok_std_h.txt
//
template <typename Container>
void stringtok (Container &container,
                std::string const &in,
                const char * const delimiters = " \t\n")
{
    const std::string::size_type len = in.length();
          std::string::size_type i = 0;

    while ( i < len )
    {
        // eat leading whitespace
        i = in.find_first_not_of (delimiters, i);

        if (i == std::string::npos)
            return;   // nothing left but white space

        // find the end of the token
        std::string::size_type j = in.find_first_of (delimiters, i);

        // push token
        if (j == std::string::npos)
        {
            container.push_back (in.substr(i));
            return;
        }
        else
            container.push_back (in.substr(i, j-i));

        // set up for next loop
        i = j + 1;
    }

    return;
}


//----------------------------------------------------------------------------------
//
// constructor: set a style during construction
//
Gnuplot::Gnuplot(bool checkForX11)
{
	this->init(checkForX11);
}


//----------------------------------------------------------------------------------
//
// constructor: set a style during construction
//
Gnuplot::Gnuplot(const std::string &style, bool checkForX11)
{
	this->init(checkForX11);
	this->set_style(style);
}


//----------------------------------------------------------------------------------
//
// constructor: open a new session, plot a signal (x)
//
Gnuplot::Gnuplot(const std::vector<double> &x,
                 const std::string &title,
                 const std::string &style,
                 const std::string &labelx,
                 const std::string &labely, bool checkForX11)
{
	this->init(checkForX11);

    this->set_style(style);
    this->set_xlabel(labelx);
    this->set_ylabel(labely);

    this->plot_x(x,title);
}


//----------------------------------------------------------------------------------
//
// constructor: open a new session, plot a signal (x,y)
//
Gnuplot::Gnuplot(const std::vector<double> &x,
                 const std::vector<double> &y,
                 const std::string &title,
                 const std::string &style,
                 const std::string &labelx,
                 const std::string &labely, bool checkForX11)
{
	this->init(checkForX11);

    this->set_style(style);
    this->set_xlabel(labelx);
    this->set_ylabel(labely);

    this->plot_xy(x,y,title);
}


//----------------------------------------------------------------------------------
//
// constructor: open a new session, plot a signal (x,y,z)
//
Gnuplot::Gnuplot(const std::vector<double> &x,
                 const std::vector<double> &y,
                 const std::vector<double> &z,
                 const std::string &title,
                 const std::string &style,
                 const std::string &labelx,
                 const std::string &labely,
                 const std::string &labelz, bool checkForX11)
{
	this->init(checkForX11);

    this->set_style(style);
    this->set_xlabel(labelx);
    this->set_ylabel(labely);
    this->set_zlabel(labelz);

    this->plot_xyz(x,y,z,title);
}


//----------------------------------------------------------------------------------
//
// Destructor: needed to delete temporary files
//
Gnuplot::~Gnuplot()
{
    if ((this->tmpfile_list).size() > 0)
    {
        for (unsigned int i = 0; i < this->tmpfile_list.size(); i++)
            remove( this->tmpfile_list[i].c_str() );

        Gnuplot::tmpfile_num -= this->tmpfile_list.size();
    }

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
    if (_pclose(this->gnucmd) == -1)
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    if (pclose(this->gnucmd) == -1)
#endif
        throw GnuplotException("Problem closing communication to gnuplot");
}


//----------------------------------------------------------------------------------
//
// Find out if valid is true
//
bool Gnuplot::is_valid()
{
    return(this->valid);
}



//----------------------------------------------------------------------------------
//
// Resets a gnuplot session (next plot will erase previous ones)
//
Gnuplot& Gnuplot::reset_plot()
{
    if (this->tmpfile_list.size() > 0)
    {
        for (unsigned int i = 0; i < this->tmpfile_list.size(); i++)
            remove(this->tmpfile_list[i].c_str());

        Gnuplot::tmpfile_num -= this->tmpfile_list.size();
        this->tmpfile_list.clear();
    }


    this->nplots = 0;
    return *this;
}


//----------------------------------------------------------------------------------
//
// Change the plotting style of a gnuplot session
//
Gnuplot& Gnuplot::set_style(const std::string &stylestr)
{
    if (stylestr.find("lines")          == std::string::npos  &&
        stylestr.find("points")         == std::string::npos  &&
        stylestr.find("linespoints")    == std::string::npos  &&
        stylestr.find("impulses")       == std::string::npos  &&
        stylestr.find("dots")           == std::string::npos  &&
        stylestr.find("steps")          == std::string::npos  &&
        stylestr.find("fsteps")         == std::string::npos  &&
        stylestr.find("histeps")        == std::string::npos  &&
        stylestr.find("labels")         == std::string::npos  &&
        stylestr.find("boxes")          == std::string::npos  &&
        stylestr.find("errorbars")      == std::string::npos  &&
        stylestr.find("xerrorbars")     == std::string::npos  &&
        stylestr.find("yerrorbars")     == std::string::npos  &&
        stylestr.find("xyerrorbars")    == std::string::npos  &&
        stylestr.find("boxerrorbars")   == std::string::npos  &&
        stylestr.find("errorlines")     == std::string::npos  &&
        stylestr.find("xerrorlines")    == std::string::npos  &&
        stylestr.find("yerrorlines")    == std::string::npos  &&
        stylestr.find("xyerrorlines")   == std::string::npos  &&
        stylestr.find("boxxyerrorbars") == std::string::npos  &&
        stylestr.find("histograms")     == std::string::npos  &&
        stylestr.find("filledcurves")   == std::string::npos  &&
        stylestr.find("financebars")    == std::string::npos  &&
        stylestr.find("candlesticks")   == std::string::npos  &&
        stylestr.find("vectors")        == std::string::npos  &&
        stylestr.find("image")          == std::string::npos  &&
        stylestr.find("rgbimage")       == std::string::npos  &&
        stylestr.find("pm3d")           == std::string::npos  )

    {
        this->pstyle = std::string("points");
    }
    else
    {
        this->pstyle = stylestr;
    }

    return *this;
}


//----------------------------------------------------------------------------------
//
// sets terminal type to windows / x11
//
Gnuplot& Gnuplot::showonscreen()
{
  this->cmd("set terminal pop");
  return *this;
}

//----------------------------------------------------------------------------------
//
// saves a gnuplot session to a postscript file
//
Gnuplot& Gnuplot::savetops(const std::string &filename)
{
  this->cmd("set terminal push");
  this->cmd("set terminal postscript color");

  std::ostringstream cmdstr;
  cmdstr << "set output \"" << filename << ".ps\"";
  this->cmd(cmdstr.str());

  return *this;
}

//----------------------------------------------------------------------------------
//
// Switches legend on
//
Gnuplot& Gnuplot::set_legend(const std::string &position)
{
    std::ostringstream cmdstr;
    cmdstr << "set key " << position;

    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
//
// Switches legend off
//
Gnuplot& Gnuplot::unset_legend()
{
    this->cmd("unset key");

    return *this;
}

//----------------------------------------------------------------------------------
//
// Turns grid on
//
Gnuplot& Gnuplot::set_grid()
{
    this->cmd("set grid");

    return *this;
}

//----------------------------------------------------------------------------------
//
// Turns grid off
//
Gnuplot& Gnuplot::unset_grid()
{
    this->cmd("unset grid");

    return *this;
}

//----------------------------------------------------------------------------------
//
// turns on log scaling for the x axis
//
Gnuplot& Gnuplot::set_xlogscale(const double base)
{
    std::ostringstream cmdstr;

    cmdstr << "set logscale x " << base;
    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
//
// turns on log scaling for the y axis
//
Gnuplot& Gnuplot::set_ylogscale(const double base)
{
    std::ostringstream cmdstr;

    cmdstr << "set logscale y " << base;
    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
//
// turns on log scaling for the z axis
//
Gnuplot& Gnuplot::set_zlogscale(const double base)
{
    std::ostringstream cmdstr;

    cmdstr << "set logscale z " << base;
    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
//
// turns off log scaling for the x axis
//
Gnuplot& Gnuplot::unset_xlogscale()
{
    this->cmd("unset logscale x");
    return *this;
}

//----------------------------------------------------------------------------------
//
// turns off log scaling for the y axis
//
Gnuplot& Gnuplot::unset_ylogscale()
{
    this->cmd("unset logscale y");
    return *this;
}

//----------------------------------------------------------------------------------
//
// turns off log scaling for the z axis
//
Gnuplot& Gnuplot::unset_zlogscale()
{
    this->cmd("unset logscale z");
    return *this;
}


//----------------------------------------------------------------------------------
//
// scales the size of the points used in plots
//
Gnuplot& Gnuplot::set_pointsize(const double pointsize)
{
    std::ostringstream cmdstr;
    cmdstr << "set pointsize " << pointsize;
    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
//
// set isoline density (grid) for plotting functions as surfaces
//
Gnuplot& Gnuplot::set_samples(const int samples)
{
    std::ostringstream cmdstr;
    cmdstr << "set samples " << samples;
    this->cmd(cmdstr.str());

    return *this;
}


//----------------------------------------------------------------------------------
//
// set isoline density (grid) for plotting functions as surfaces
//
Gnuplot& Gnuplot::set_isosamples(const int isolines)
{
    std::ostringstream cmdstr;
    cmdstr << "set isosamples " << isolines;
    this->cmd(cmdstr.str());

    return *this;
}


//----------------------------------------------------------------------------------
//
// enables hidden line removal for surface plotting
//
Gnuplot& Gnuplot::set_hidden3d()
{
    this->cmd("set hidden3d");

    return *this;
}

//----------------------------------------------------------------------------------
//
// disables hidden line removal for surface plotting
//
Gnuplot& Gnuplot::unset_hidden3d()
{
    this->cmd("unset hidden3d");

    return *this;
}

//----------------------------------------------------------------------------------
//
// Sets the title of a gnuplot session
//
Gnuplot& Gnuplot::set_title(const std::string &title)
{
    std::ostringstream cmdstr;

    cmdstr << "set title \"" << title << "\"";
    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
//
// Clears the title of a gnuplot session
//
Gnuplot& Gnuplot::unset_title()
{
    this->set_title("");

    return *this;
}

//----------------------------------------------------------------------------------
//
// set labels
//
// set the xlabel
Gnuplot& Gnuplot::set_xlabel(const std::string &label)
{
    std::ostringstream cmdstr;

    cmdstr << "set xlabel \"" << label << "\"";
    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
// set the ylabel
//
Gnuplot& Gnuplot::set_ylabel(const std::string &label)
{
    std::ostringstream cmdstr;

    cmdstr << "set ylabel \"" << label << "\"";
    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
// set the zlabel
//
Gnuplot& Gnuplot::set_zlabel(const std::string &label)
{
    std::ostringstream cmdstr;

    cmdstr << "set zlabel \"" << label << "\"";
    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
//
// set range
//
// set the xrange
Gnuplot& Gnuplot::set_xrange(const int iFrom,
                             const int iTo)
{
    std::ostringstream cmdstr;

    cmdstr << "set xrange[" << iFrom << ":" << iTo << "]";
    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
// set the yrange
//
Gnuplot& Gnuplot::set_yrange(const int iFrom,
                             const int iTo)
{
    std::ostringstream cmdstr;

    cmdstr << "set yrange[" << iFrom << ":" << iTo << "]";
    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
// set the zrange
//
Gnuplot& Gnuplot::set_zrange(const int iFrom,
                             const int iTo)
{
    std::ostringstream cmdstr;

    cmdstr << "set zrange[" << iFrom << ":" << iTo << "]";
    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
//
// set the palette range
//
Gnuplot& Gnuplot::set_cbrange(const int iFrom,
                              const int iTo)
{
    std::ostringstream cmdstr;

    cmdstr << "set cbrange[" << iFrom << ":" << iTo << "]";
    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
//
// Plots a linear equation y=ax+b (where you supply the
// slope a and intercept b)
//
Gnuplot& Gnuplot::plot_slope(const double a,
                             const double b,
                             const std::string &title)
{
    std::ostringstream titlestr;
    std::ostringstream cmdstr;

    if (title == "")
        titlestr << "f(x) = " << a << " * x + " << b;
    else
        titlestr << title;

    if (this->nplots > 0  &&  this->two_dim == true)
        cmdstr << "replot " << a << " * x + " << b << " title \""
               << titlestr.str() << "\" with " << this->pstyle;
    else
        cmdstr << "plot " << a << " * x + " << b << " title \""
               << titlestr.str() << "\" with " << this->pstyle;

    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
//
// Plot an equation which is supplied as a std::string y=f(x) (only f(x) expected)
//
Gnuplot& Gnuplot::plot_equation(const std::string &equation,
                                const std::string &title)
{
    std::string titlestr, plotstr;
    std::ostringstream cmdstr;

    if (title == "")
        titlestr = "f(x) = " + equation;
    else
        titlestr = title;

    if (this->nplots > 0  &&  this->two_dim == true)
        plotstr = "replot";
    else
        plotstr = "plot";

    cmdstr << plotstr << " " << equation << " " << "title \""
           << titlestr << "\" with " << this->pstyle;

    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
//
// plot an equation supplied as a std::string y=(x)
//
Gnuplot& Gnuplot::plot_equation3d(const std::string &equation,
                                  const std::string &title)
{
    std::string titlestr, plotstr;
    std::ostringstream cmdstr;

    if (title == "")
        titlestr = "f(x,y) = " + equation;
    else
        titlestr = title;

    if (this->nplots > 0  &&  this->two_dim == false)
        plotstr = "replot";
    else
        plotstr = "splot";

    cmdstr << plotstr << " " << equation << " " << "title \""
           << titlestr << "\" with " << this->pstyle;

    this->cmd(cmdstr.str());

    return *this;
}


//----------------------------------------------------------------------------------
//
// Plots a 2d graph from a list of doubles (x) saved in a file
//
Gnuplot& Gnuplot::plotfile_x(const std::string &filename,
                             const std::string &title)
{
    //
    // check if file exists
    //
    if( !(Gnuplot::file_exists(filename,4)) ) // check existence and read permission
    {
        std::ostringstream except;
        except << "File \"" << filename << "\" does not exist";
        throw GnuplotException( except.str() );
        return *this;
    }


    std::ostringstream cmdstr;
    //
    // command to be sent to gnuplot
    //
    if (this->nplots > 0  &&  this->two_dim == true)
        cmdstr << "replot ";
    else
        cmdstr << "plot ";

    if (title == "")
        cmdstr << "\"" << filename << "\" notitle with " << this->pstyle;
    else
        cmdstr << "\"" << filename << "\" title \"" << title << "\" with " << this->pstyle;

    //
    // Do the actual plot
    //
    this->cmd(cmdstr.str()); //nplots++; two_dim = true;  already in this->cmd();

    return *this;
}


//----------------------------------------------------------------------------------
//
// Plots a 2d graph from a list of doubles: x
//
Gnuplot& Gnuplot::plot_x(const std::vector<double> &x,
                         const std::string &title)
{
    if (x.size() == 0)
    {
        throw GnuplotException("std::vector too small");
        return *this;
    }

    std::ofstream tmp;
    std::string name = create_tmpfile(tmp);
    if (name == "")
        return *this;

    //
    // write the data to file
    //
    for (unsigned int i = 0; i < x.size(); i++)
        tmp << x[i] << std::endl;

    tmp.flush();
    tmp.close();


    this->plotfile_x(name, title);

    return *this;
}


//----------------------------------------------------------------------------------
//
// Plots a 2d graph from a list of doubles (x y) saved in a file
//
Gnuplot& Gnuplot::plotfile_xy(const std::string &filename,
                              const std::string &title)
{
    //
    // check if file exists
    //
    if( !(Gnuplot::file_exists(filename,4)) ) // check existence and read permission
    {
        std::ostringstream except;
        except << "File \"" << filename << "\" does not exist";
        throw GnuplotException( except.str() );
        return *this;
    }


    std::ostringstream cmdstr;
    //
    // command to be sent to gnuplot
    //
    if (this->nplots > 0  &&  this->two_dim == true)
        cmdstr << "replot ";
    else
        cmdstr << "plot ";

    if (title == "")
        cmdstr << "\"" << filename << "\" notitle with " << this->pstyle;
    else
        cmdstr << "\"" << filename << "\" title \"" << title << "\" with " << this->pstyle;

    //
    // Do the actual plot
    //
    this->cmd(cmdstr.str());

    return *this;
}

//----------------------------------------------------------------------------------
//
// Plots a 2d graph from a list of doubles: x y
//
Gnuplot& Gnuplot::plot_xy(const std::vector<double> &x,
                          const std::vector<double> &y,
                          const std::string &title)
{
    if (x.size() == 0 || y.size() == 0)
    {
        throw GnuplotException("std::vectors too small");
        return *this;
    }

    if (x.size() != y.size())
    {
        throw GnuplotException("Length of the std::vectors differs");
        return *this;
    }


    std::ofstream tmp;
    std::string name = create_tmpfile(tmp);
    if (name == "")
        return *this;

    //
    // write the data to file
    //
    for (unsigned int i = 0; i < x.size(); i++)
        tmp << x[i] << " " << y[i] << std::endl;

    tmp.flush();
    tmp.close();


    this->plotfile_xy(name, title);

    return *this;
}


//----------------------------------------------------------------------------------
//
// Plots a 2d graph with errorbars from a list of doubles (x y dy) saved in a file
//
Gnuplot& Gnuplot::plotfile_xy_err(const std::string &filename,
                                  const std::string &title)
{
    //
    // check if file exists
    //
    if( !(Gnuplot::file_exists(filename,4)) ) // check existence and read permission
    {
        std::ostringstream except;
        except << "File \"" << filename << "\" does not exist";
        throw GnuplotException( except.str() );
        return *this;
    }


    std::ostringstream cmdstr;
    //
    // command to be sent to gnuplot
    //
    if (this->nplots > 0  &&  this->two_dim == true)
        cmdstr << "replot ";
    else
        cmdstr << "plot ";

    if (title == "")
    {
        if (this->pstyle.compare(0, 9, "errorbars") != 0)
            cmdstr << "\"" << filename << "\" notitle with " << this->pstyle << ", \""
                   << filename << "\" notitle with errorbars";
        else
            cmdstr << "\"" << filename << "\" notitle with errorbars";

    }
    else
    {
        if (this->pstyle.compare(0, 9, "errorbars") != 0)
            cmdstr << "\"" << filename << "\" title \"" << title << "\" with "
                   << this->pstyle << ", \"" << filename << "\" notitle with errorbars";

        else
            cmdstr << "\"" << filename << "\" title \"" << title << "\" with errorbars";
    }

    //
    // Do the actual plot
    //
    this->cmd(cmdstr.str());

    return *this;
}


//----------------------------------------------------------------------------------
//
// plot x,y pairs with dy errorbars
//
Gnuplot& Gnuplot::plot_xy_err(const std::vector<double> &x,
                              const std::vector<double> &y,
                              const std::vector<double> &dy,
                              const std::string &title)
{
    if (x.size() == 0 || y.size() == 0 || dy.size() == 0)
    {
        throw GnuplotException("std::vectors too small");
        return *this;
    }

    if (x.size() != y.size() || y.size() != dy.size())
    {
        throw GnuplotException("Length of the std::vectors differs");
        return *this;
    }


    std::ofstream tmp;
    std::string name = create_tmpfile(tmp);
    if (name == "")
        return *this;

    //
    // write the data to file
    //
    for (unsigned int i = 0; i < x.size(); i++)
        tmp << x[i] << " " << y[i] << " " << dy[i] << std::endl;

    tmp.flush();
    tmp.close();


    // Do the actual plot
    this->plotfile_xy_err(name, title);

    return *this;
}


//----------------------------------------------------------------------------------
//
// Plots a 3d graph from a list of doubles (x y z) saved in a file
//
Gnuplot& Gnuplot::plotfile_xyz(const std::string &filename,
                               const std::string &title)
{
    //
    // check if file exists
    //
    if( !(Gnuplot::file_exists(filename,4)) ) // check existence and read permission
    {
        std::ostringstream except;
        except << "File \"" << filename << "\" does not exist";
        throw GnuplotException( except.str() );
        return *this;
    }


    std::ostringstream cmdstr;
    //
    // command to be sent to gnuplot
    //
    if (this->nplots > 0  &&  this->two_dim == false)
        cmdstr << "replot ";
    else
        cmdstr << "splot ";

    if (title == "")
        cmdstr << "\"" << filename << "\" notitle with " << this->pstyle;
    else
        cmdstr << "\"" << filename << "\" title \"" << title << "\" with " << this->pstyle;

    //
    // Do the actual plot
    //
    this->cmd(cmdstr.str());

    return *this;
}


//----------------------------------------------------------------------------------
//
// Plots a 3d graph from a list of doubles: x y z
//
Gnuplot& Gnuplot::plot_xyz(const std::vector<double> &x,
                           const std::vector<double> &y,
                           const std::vector<double> &z,
                           const std::string &title)
{
    if (x.size() == 0 || y.size() == 0 || z.size() == 0)
    {
        throw GnuplotException("std::vectors too small");
        return *this;
    }

    if (x.size() != y.size() || x.size() != z.size())
    {
        throw GnuplotException("Length of the std::vectors differs");
        return *this;
    }


    std::ofstream tmp;
    std::string name = create_tmpfile(tmp);
    if (name == "")
        return *this;

    //
    // write the data to file
    //
    for (unsigned int i = 0; i < x.size(); i++)
    {
        tmp << x[i] << " " << y[i] << " " << z[i] <<std::endl;
    }

    tmp.flush();
    tmp.close();


    this->plotfile_xyz(name, title);

    return *this;
}



//----------------------------------------------------------------------------------
//
/// *  note that this function is not valid for versions of GNUPlot below 4.2
//
Gnuplot& Gnuplot::plot_image(const unsigned char * ucPicBuf,
                             const int iWidth,
                             const int iHeight,
                             const std::string &title)
{
    std::ofstream tmp;
    std::string name = create_tmpfile(tmp);
    if (name == "")
        return *this;

    //
    // write the data to file
    //
    int iIndex = 0;
    for(int iRow = 0; iRow < iHeight; iRow++)
    {
        for(int iColumn = 0; iColumn < iWidth; iColumn++)
		{
			tmp << iColumn << " " << iRow  << " " << static_cast<float>(ucPicBuf[iIndex++]) << std::endl;
	  	}
    }

    tmp.flush();
    tmp.close();


    std::ostringstream cmdstr;
    //
    // command to be sent to gnuplot
    //
    if (this->nplots > 0  &&  this->two_dim == true)
        cmdstr << "replot ";
    else
        cmdstr << "plot ";

    if (title == "")
        cmdstr << "\"" << name << "\" with image";
    else
        cmdstr << "\"" << name << "\" title \"" << title << "\" with image";

    //
    // Do the actual plot
    //
    this->cmd(cmdstr.str());

    return *this;
}



//----------------------------------------------------------------------------------
//
// Sends a command to an active gnuplot session
//
Gnuplot& Gnuplot::cmd(const std::string &cmdstr)
{
    if( !(this->valid) )
    {
        return *this;
    }

    //std::string tmp = cmdstr + "\n";

    // writes the string cmdstr to the stream gnucmd
    fputs( (cmdstr+"\n").c_str(), this->gnucmd );

    // flushes the stream, the stream remains open after this call
    fflush(this->gnucmd);


    if( cmdstr.find("replot") != std::string::npos )
    {
        this->nplots++;
    }
    else if( cmdstr.find("splot") != std::string::npos )
    {
        this->two_dim = false;
        this->nplots++;
    }
    else if( cmdstr.find("plot") != std::string::npos )
    {
        this->two_dim = true;
        this->nplots++;
    }

    return *this;
}



//----------------------------------------------------------------------------------
//
// Sends a command to an active gnuplot session, identical to cmd()
//
Gnuplot& Gnuplot::operator<<(const std::string &cmdstr)
{
    this->cmd(cmdstr);
    return *this;
}



//----------------------------------------------------------------------------------
//
// Opens up a gnuplot session, ready to receive commands
//
void Gnuplot::init(bool checkForX11)
{
  bool didCheck;
  didCheck=false;
#if defined(unix) || defined(__unix) || defined(__unix__)  // check DISPLAY variable
  if (getenv("DISPLAY") == NULL)   {
    this->valid = false;
    throw GnuplotException("Can't find DISPLAY variable");
  }
#elif defined(__APPLE__)
  if (checkForX11) {
    if (getenv("DISPLAY") == NULL)   {
      this->valid = false;
      throw GnuplotException("Can't find DISPLAY variable");
    }   
    didCheck = true;
  } else {
    this->valid = true;
    didCheck = false;
  }
#endif

    // if gnuplot not available
    if (!Gnuplot::get_program_path())
    {
        this->valid = false;
        throw GnuplotException("Can't find gnuplot");
    }

    //
    // open pipe
    //
    std::string tmp = Gnuplot::m_sGNUPlotPath + "/" + m_sGNUPlotFileName;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
    this->gnucmd = _popen(tmp.c_str(),"w");
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    this->gnucmd = popen(tmp.c_str(),"w");
#endif

    if (!this->gnucmd)
    {
        this->valid = false;
        throw GnuplotException("Couldn't open connection to gnuplot");
    }

    this->nplots = 0;
    this->valid = true;
    
    // set terminal type
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
    this->cmd("set terminal windows");
#elif defined(unix) || defined(__unix) || defined(__unix__)
    this->cmd("set terminal x11");
#elif defined(__APPLE__)
    if (didCheck) this->cmd("set terminal x11");
    else this->cmd("set terminal aqua");
#endif


    return;
}


//----------------------------------------------------------------------------------
//
// Find out if a command lives in m_sGNUPlotPath or in PATH
//
bool Gnuplot::get_program_path()
{
    //
    // first look in m_sGNUPlotPath for Gnuplot
    //
    std::string tmp = Gnuplot::m_sGNUPlotPath + "/" + Gnuplot::m_sGNUPlotFileName;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
    if ( Gnuplot::file_exists(tmp,0) ) // check existence
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    if ( Gnuplot::file_exists(tmp,1) ) // check existence and execution permission
#endif
    {
        return true;
    }


    //
    // second look in PATH for Gnuplot
    //
    char *path;
    // Retrieves a C string containing the value of the environment variable PATH
    path = getenv("PATH");


    if (path == NULL)
    {
        throw GnuplotException("Path is not set");
        return false;
    }
    else
    {
        std::list<std::string> ls;

        //split path into list ls
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
        stringtok(ls,path,";");
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
        stringtok(ls,path,":");
#endif

        // scan list for Gnuplot program files
        for (std::list<std::string>::const_iterator i = ls.begin(); i != ls.end(); ++i)
        {
            tmp = (*i) + "/" + Gnuplot::m_sGNUPlotFileName;
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
            if ( Gnuplot::file_exists(tmp,0) ) // check existence
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
            if ( Gnuplot::file_exists(tmp,1) ) // check existence and execution permission
#endif
            {
                Gnuplot::m_sGNUPlotPath = *i; // set m_sGNUPlotPath
                return true;
            }
        }

        tmp = "Can't find gnuplot neither in PATH nor in \"" + Gnuplot::m_sGNUPlotPath + "\"";
        throw GnuplotException(tmp);

        Gnuplot::m_sGNUPlotPath = "";
        return false;
    }
}



//----------------------------------------------------------------------------------
//
// check if file exists
//
bool Gnuplot::file_exists(const std::string &filename, int mode)
{
    if ( mode < 0 || mode > 7)
    {
        throw std::runtime_error("In function \"Gnuplot::file_exists\": mode has to be an integer between 0 and 7");
        return false;
    }

    // _access(const char *path, int mode) returns 0 if the file has the given mode,
    // it returns -1 if the named file does not exist or is not accessible in the given mode
    // mode = 0 (F_OK) (default): checks file for existence only
    // mode = 1 (X_OK): execution permission
    // mode = 2 (W_OK): write permission
    // mode = 4 (R_OK): read permission
    // mode = 6       : read and write permission
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
    if (_access(filename.c_str(), mode) == 0)
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    if (access(filename.c_str(), mode) == 0)
#endif
    {
        return true;
    }
    else
    {
        return false;
    }

}



//----------------------------------------------------------------------------------
//
// Opens a temporary file
//
std::string Gnuplot::create_tmpfile(std::ofstream &tmp)
{

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
    char name[] = "gnuplotiXXXXXX"; //tmp file in working directory
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    char name[] = "/tmp/gnuplotiXXXXXX"; // tmp file in /tmp
#endif

    //
    // check if maximum number of temporary files reached
    //
    if (Gnuplot::tmpfile_num == GP_MAX_TMP_FILES - 1)
    {
        std::ostringstream except;
        except << "Maximum number of temporary files reached (" << GP_MAX_TMP_FILES
               << "): cannot open more files" << std::endl;

        throw GnuplotException( except.str() );
        return "";
    }

    //
    // open temporary files for output
    //
#if defined(_MSC_VER) && ( defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__) )
    if (_mktemp_s(name) == NULL)
#elif defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
    if (_mktemp(name) == NULL)
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    if (mkstemp(name) == -1)
#endif
    {
        std::ostringstream except;
        except << "Cannot create temporary file \"" << name << "\"";
        throw GnuplotException(except.str());
        return "";
    }

    tmp.open(name);
    if (tmp.bad())
    {
        std::ostringstream except;
        except << "Cannot create temporary file \"" << name << "\"";
        throw GnuplotException(except.str());
        return "";
    }

    //
    // Save the temporary filename
    //
    this->tmpfile_list.push_back(name);
    Gnuplot::tmpfile_num++;

    return name;
}
