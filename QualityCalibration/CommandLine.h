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

#ifndef COMMANDLINE_H_
#define COMMANDLINE_H_

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <sstream>
#include <stdexcept>
#include "stringify.h"

using namespace std;

class CommandLine {
public:

  typedef map<string,string> parm_types_t;
  typedef map<string,string> parm_description_t;
  typedef map<string,bool>   parm_required_t;

private:
	
  // Constructor is private. It can only be called from the 'Instance' class method.
  // That method only calls it once, to create a single instance, and returns a pointer
  // to the single instance on subsequent calls.
  
  CommandLine (ostream &err_in=std::cerr) : err(err_in) {
  }
	
public:

  // This is the class method a user calls to get a pointer to the unique instance
  // of a CommandLine object.
  
  static CommandLine *Instance () {
    if (only_instance == 0) {
      only_instance = new CommandLine;
    }
    return only_instance;
  }
  
  // Set valid parameters, with descriptions and default values. The next steps 
  // are expected to be adding command-line parameters, and adding parameters from
  // a config file. Command line parameters override those specified in the config
  // file. But you need to process the command line to find the name of the config
  // file. That means default parameters have to be saved in a separate map, and
  // applied at the end (in CommandLine::check), rather than the more straightforward
  // approach of stuffing them into parm_list and then overwriting them.
  
  void add_valid_parm (const string &parm,
                       const string &description,
                       bool required=false,
                       string default_value="") {
    
    parm_description[parm] = description;
    parm_required[parm]    = required;
    
    if (default_value != "") {
      if (required) {
        throw (std::logic_error("required parameters don't need defaults: " + parm));
      }
      parm_default[parm] = default_value;
    }
  
  }

  // Parse command line parameters of the form "--parm value" (where value is optional) 
  // into the parm_list map. No parms_types checking is done here -- that's in check(),
  // so it can be invoked separately. Stop when the next parameter doesn't begin with '--'.
  // Return the index of the first parameter in argv which wasn't processed.
  
  // There's a semantic difficulty here with parameters which are true/false flags.
  // A parameter specified with no value is interpreted as having the value "1" (which
  // converts to true when treated as a bool). The string "true" and "yes" will be saved
  // as "1", "false" and "no" as "0". This allows a parameter to be defaulted true and 
  // turned off (false) -- either in the config file, or by saying "--fizzle false" on
  // the command line.

  // TODO: Accept parameters of the form "--parm=value"

	int process_args (int argc, char **argv) {

		vector<string> argvec; // string version of argv

		// index from zero to preserve original numbering
		for (int i=0; i<argc; ++i) {
			argvec.push_back(string(argv[i]));
		}

		int i = 1;		// index from 1 to skip program name
		while (i<argc && argvec[i].substr(0, 2) == "--") {
		  
		  string parm (argvec[i].substr(2));
		  
		  // Check whether parm is valid (i.e., was specified via add_valid_parm). We'll use
		  // parm_description for that, as it's updated by every call to add_valid_parm.
		  
		  if (parm_description.count(parm) == 0) {
        throw (std::invalid_argument("invalid command line parameter: " + parm));
		  }
		  
		  string value;
		  
			if (i+1 < argc && argvec[i+1].substr(0, 2) != "--") {    // is next token a value?
				value = argvec[i+1];
				i += 2;         // eat parm and value
			} else {
				value = "1";
				++i;            // eat parm only
			}

      if (value == "true" || value == "yes") {
        value = "1";
      } else if (value == "false" || value == "no") {
        value = "0";
      }

      parm_list[parm] = value;
			
		}

		return i;

	}

	// Read parameters from a specified file and add them to parm_list as if they were 
	// read from the command line. Config file is expected to contain lines containing
	// two fields, a parameter name and a value, separated by spaces. Fields beyond the
	// first two on a line are ignored. Lines beginning with # considered comments and
	// ignored. Blank lines are ignored.

	bool read_config_file (const string &filename) {

	  ifstream in(filename.c_str());
	  string line;
	  istringstream fields;
	  string parm;
	  string value;
	  
	  while (getline (in, line)) {

	    if (line.find_first_not_of(' ') != string::npos
	        && line.substr(0, 1) != "#") {
	    
	      fields.str(line);
	      fields >> parm;
	      fields >> value;
	      fields.clear();     // clear eof flag
	      
	      // See comment block in process_args
	      if (value.empty() || value == "true" || value == "yes") {
	        value = "1";
	      } else if (value == "false" || value == "no") {
	        value = "0";
	      }

	      err << "parm: " << parm << " value: |" << value << "|" << endl; 
	      
	      if (parm_description.count(parm) == 0) {
	        throw (std::invalid_argument("invalid config file parameter: " + parm));
	      }
	    
	      if (parm_list.count(parm) == 0) {    // if not already there from command line 
	        parm_list[parm] = value;
	      }
	      
	    }

	  }
	  
	  return true;
	  
	}

	// Add default values, check that required parameters are present.
	
	bool check () {

	  // Stuff in default values for parameters that haven't yet been specified.
	  
    for(parm_types_t::const_iterator i=parm_default.begin(); i != parm_default.end(); i++) {
      
      if (parm_list.count((*i).first) == 0) {

        // For consistency's sake, we'll do the true/false/yes/no conversion to boolean
        // (as described in the comment block above process_args) for the default values as 
        // well. We will *not* convert an empty string to "1" (as is done for command line
        // parameters), however, since the caller probably really wanted an empty string.
        
        string value ((*i).second);
        
        if (value == "true" || value == "yes") {
          value = "1";
        } else if (value == "false" || value == "no") {
          value = "0";
        }

        parm_list[(*i).first] = value;

      }
    }

    // Check for presence of required parameters.
    
    for(parm_required_t::const_iterator i=parm_required.begin(); i != parm_required.end(); i++) {
      if ((*i).second) {
        if (parm_list.count((*i).first) == 0) {
          throw (std::invalid_argument("missing required parameter: " + (*i).first));
        } 
      } 
    }
	  
    return true;
	  
	}

	// Return true if parameter is set. Whether it was set by command line, config file,
	// or default (an odd thing to do) makes no difference.

	bool is_set (const string &parm) {
		return (parm_list.count(parm) > 0);
	}

	string get_parm (const string &parm) {
    // no bounds checking...
		if(is_set(parm)) {
			return (*parm_list.find(parm)).second;
		} else {
			return string("");
		}
	}

	template <class _T>
	_T get_parm_as (const string &parm) {
	  // no bounds checking...
	  return convertTo<_T>(get_parm(parm));
	}

	/// Print usage information
	string usage() {

		ostringstream ss;
		for(parm_description_t::const_iterator i=parm_description.begin();i != parm_description.end();i++) {
			ss << setw(40) << (*i).first << " : " << (*i).second << endl;
		}

		return ss.str();
	}

  string dump_settings() {
    ostringstream ss;
    for(parm_description_t::const_iterator i=parm_description.begin();i != parm_description.end();i++) {
      ss << setw(40) << (*i).first << " : " << get_parm((*i).first) << endl;
    }

    return ss.str();      
  }
  
  string dump_settings_xml() {
    ostringstream ss;

    ss << "<parameters>" << endl;
    for(parm_description_t::const_iterator i=parm_description.begin();i != parm_description.end();i++) {
      ss << "<parameter name=\"" << (*i).first << "\" value=\"" <<  get_parm((*i).first) << "\" />" << endl;
    }
    ss << "</parameters>" << endl;

    return ss.str();      
  }

private:
  
	parm_types_t         parm_list;        ///< Parameter values
	parm_description_t   parm_description; ///< Description of valid parameters/usage
	parm_required_t      parm_required;    ///< flags indicating whether parm is required
	parm_types_t         parm_default;     ///< default values for parms, specified using add_valid_parm
	ostream              &err;             ///< Error output will be writen here, set to cerr in constructor default
	
	static CommandLine   *only_instance;
	
};

CommandLine   *CommandLine::only_instance;


#endif /*COMMANDLINE_H_*/
