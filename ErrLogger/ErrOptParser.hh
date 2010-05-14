#ifndef ERROPTPARSER_HH
#define ERROPTPARSER_HH

//--------------------------------------------------------------------------
//
// $Id: ErrOptParser.hh,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
//
// Environment:
//      This software was developed for the BaBar collaboration.  If you
//      use all or part of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 2002   California Institute of Technology
//
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------
#include <string>
#include <map>
#include <iostream>

#include "ErrLogger/ErrOpt.hh"

/**
 *  ErrOptParser keeps a register of acceptable error suboptions for a given
 *  derivative of ErrLog and processes values specified at the command 
 *  line.
 *
 *  How to use it:
 *    1) Create an ErrOptParser.
 *    2) Instantiate some variables to hold suboption values.
 *    3) Instantiate ErrOptRef<T>'s specifying names and value-parsing
 *       behavior for suboptions and attach them to the corresponding
 *       variables.
 *    4) Call registerOption() and registerNegatableOption() as desired.
 *    5) Obtain a suboption string suitable for use as the first 
 *       argument to Unix98 getsubopt().  Typically this might be the
 *       value for the "-E" option on the command line.
 *    6) Pass this string to parseOpt() and check return value for success.
 *    7) Use the values in the option variables.
 *
 *  The class is intended to be extensible so that additional types of
 *  options beyond the two already provided can be defined.  (This hasn't
 *  been attempted; please complain if you try to do so and the right 
 *  things haven't been exposed as "protected".
 *
 *  This software was developed for the BaBar collaboration.  If you
 *  use all or part of it, please give an appropriate acknowledgement.
 *
 *  Copyright (C) 2002   California Institute of Technology
 *
 *  @version $Id: ErrOptParser.hh,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
 *
 *  @author Matt Dorsten  (originator)
 *  @author Gregory Dubois-Felsmann
 */

class ErrOptParser {

public:
  ErrOptParser() : optmap_() { }
  virtual ~ErrOptParser(); 

  /**
   *  Register a named suboption associated with a value of type T.
   *  Instances of "<name>=<value>" in the parsed suboption string
   *  will result in the value string being converted to T using the
   *  operator>>(std::istream&, T&) stream extractor.  An error message and
   *  failure status will be returned if that conversion fails.
   *
   *  Whether the value is required to be present or not is controlled
   *  by the ErrOpt::ValType used in the construction of the ErrOptRef<T>.
   *
   *  Following parsing, the ErrOptRef<T>::present() function can be used
   *  to determine whether the suboption was specified at all.
   */
  template<class T> bool registerOption( ErrOptRef<T>& );

  /**
   *  Register a named negatable Boolean suboption.  This causes the
   *  following forms to be acceptable in the parse:
   *
   *     <name>
   *     <name>=<value>
   *     no<name>
   *
   *  Whether one or the other or both of the first two are accepted 
   *  is under the the control of the ErrOpt::ValType used in the
   *  construction of the ErrOptRef<bool>.
   *
   *  The acceptable value strings are "true", "1", "false", and "0",
   *  the ones understood when "boolalpha" is used in Standard iostreams
   *  input of bool values (i.e., not "yes" or "no" and not abbreviations),
   *  for consistency with the templated form of registerOption().  This 
   *  is subject to future revision.
   *
   *  An error message and failure status will be returned if this conversion
   *  fails.
   *
   *  Following parsing, the ErrOptRef<bool>::present() function can be 
   *  used to determine whether the suboption was specified at all.  The
   *  presence of either the "<name>" or "no<name>" form will cause
   *  present() to return true.
   */
  bool registerNegatableOption( ErrOptRef<bool>& );

  /** 
   *  Parse a suboption string suitable for use as the first argument to
   *  Unix98 getsubopt().  The string is copied internally before being
   *  passed to getsubopt(), however, so the caller's string will not see
   *  the modifications (insertions of '\0') that getsubopt() performs on
   *  its input.
   *
   *  Returns false if any unregistered option is found, if the 
   *  ErrOpt::ValType constraint for any registered option is violated,
   *  or if the parsing of any value fails.  It is unspecified whether
   *  all the parsing will have been attempted if any of it fails.  It is
   *  intended that an application not try to proceed beyond such a failure
   *  and use the variables attached to the options; they may not be in a
   *  consistent state.
   */
  virtual bool parseOpts( const char* );
  
protected:

  /**
   *  Process the value for a single named suboption.  The second 
   *  argument is intended to be the value pointer obtained from 
   *  getsubopt(), and so can be 0 if the suboption didn't come with
   *  an "=<value>" string in the input.
   */
  virtual bool processOption( const std::string&, const char* );

  class OptHelperBase{
  protected:
    OptHelperBase( ) { }
  public:
    virtual ~OptHelperBase( ) { }
    virtual bool processValue(const char*) = 0;
  };  

  template<class T> class OptHelper : public OptHelperBase{
  public:
    OptHelper( ErrOptRef<T>& optref )
      : OptHelperBase(), optref_(optref) { }
    virtual bool processValue(const char*);
  private:
    ErrOptRef<T>& optref_;
  };
  
  class BoolHelper : public OptHelperBase{
  public:
    BoolHelper( ErrOptRef<bool>& optref, bool negated )
      : OptHelperBase(), optref_(optref), negated_(negated) { }
    virtual bool processValue(const char*);
  private:
    // need to fix when Standard iostreams bool extraction is available  FIXME
    static bool getbool(const std::string&, bool&);

    ErrOptRef<bool>& optref_;
    bool negated_;
  };

  // The nested classes above need to be made friends so that they can use
  // the "setPresent" door into the private interface of ErrOpt.  Note
  // that this is not appropriate to do through the OptHelperBase class,
  // because not all possible subclasses of it might concern themselves
  // with ErrOpt objects at all.
  template <class T> friend class OptHelper;
  friend class BoolHelper;

  static void setPresent( ErrOpt& eo, bool b ) { eo.setPresent(b); }


private:
  ErrOptParser(const ErrOptParser&);            // prevent copying
  ErrOptParser& operator=(const ErrOptParser&);

  typedef std::map<std::string, OptHelperBase*> MapType;
  MapType optmap_;
};
  
#endif
