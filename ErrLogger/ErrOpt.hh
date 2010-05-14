#ifndef ERROPT_HH
#define ERROPT_HH

//--------------------------------------------------------------------------
//
// $Id: ErrOpt.hh,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
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

//		---------------------
// 		-- Class Interface --
//		---------------------

/**
 *  This class represents a single error suboption: its name, whether
 *  it has a required, optional, or no value associated with it, and
 *  whether the suboption has actually been found in parsing a suboption
 *  list.
 *
 *  The templated subclass ErrOptRef<T> attaches an ErrOpt to a reference 
 *  to a variable of type T that is intended to receive the suboption value.
 *
 *  This software was developed for the BaBar collaboration.  If you
 *  use all or part of it, please give an appropriate acknowledgement.
 *
 *  Copyright (C) 2002   California Institute of Technology
 *
 *  @see TemplateTemplate
 *
 *  @version $Id: ErrOpt.hh,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
 *
 *  @author Matt Dorsten  (originator)
 */

class ErrOpt {
public:
  /**
   *  Based on the Unix98 getsubopt() behavior for parsing suboptions,
   *  where each one is expected to have the form "<name>[=<value>]",
   *  this enum controls whether a given suboption's associated value
   *  will be treated as required, optional, or prohibited when parsing
   *  the string.  ("Prohibited" is probably only useful with negatable
   *  boolean options.)
   *
   *  @see ErrOptParser::registerNegatable
   */
  enum ValType { Required, Optional, Prohibited };

protected:
  ErrOpt( const std::string& name, ValType vt ) 
    : present_(false), name_(name), type_(vt) { }

public:
  const std::string& name() const { return name_; }
  ValType type() const { return type_; }
  bool present() const { return present_; }

private:
  // Only the parser is permitted to flag an option as present, so access
  // to this interface is limited:
  friend class ErrOptParser;
  void setPresent( bool b ) { present_ = b; }

  bool present_;

  const std::string name_;
  const ValType type_;
};


/**
 *  This class attaches an ErrOpt to a reference to a variable of type T 
 *  that is intended to receive the suboption value.
 */

template <class T>
class ErrOptRef : public ErrOpt {
public:
  ErrOptRef( const std::string& name, ValType vt, T& t )
    : ErrOpt(name, vt), t_(t) { }
  T& value() { return t_; }
private:
  T& t_;
};

#endif // ERROPT_HH
