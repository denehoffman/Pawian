//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: ErrOptParser.cc,v 1.2 2006/04/19 15:59:06 steinke Exp $
//
// Description:
//	ErrOptParser implementation
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Matt Dorsten  (originator)
//
// Copyright Information:
//	Copyright (C) 2002   California Institute of Technology
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "ErrLogger/ErrOptParser.hh"

//-----------------
// C/C++ Headers --
//-----------------
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstring>
#include <stdlib.h>
using std::cerr;
using std::endl;
using std::istringstream;
using std::stringstream;

ErrOptParser::~ErrOptParser()
{
  for(MapType::iterator opti = optmap_.begin(); opti != optmap_.end(); ++opti){
    delete opti->second;
    opti->second = 0;
  }
}

template<class T> 
bool 
ErrOptParser::registerOption( ErrOptRef<T>& optref )
{
  MapType::iterator opti = optmap_.find(optref.name());
  if(opti != optmap_.end()){
    cerr << "Duplicate suboption " << optref.name() << endl;
    return false;
  }
  else
    optmap_[optref.name()] = new OptHelper<T>(optref);
  return true;
}

// must register only the positive option, not its negation
bool 
ErrOptParser::registerNegatableOption( ErrOptRef<bool>& optref )
{
  std::string noname("no");
  noname += optref.name();

  MapType::iterator opti = optmap_.find(optref.name());
  MapType::iterator noopti = optmap_.find(noname);
  if(opti != optmap_.end()){
    cerr << "Duplicate suboption \"" << optref.name() << "\"" << endl;
    return false;
  }
  else if(noopti != optmap_.end()){
    cerr << "Duplicate suboption \"" << noname << "\"" << endl;
    return false;
  }
  else{
    optmap_[optref.name()] = new BoolHelper(optref, true);
    optmap_[noname]        = new BoolHelper(optref, false);
  }
  return true;
}

// This class is used to make a writable copy of a string that
// getsubopt can modify, and that is represented as a char*, while
// retaining ownership of it.
struct StringCopy {
  char* ptr_;
  explicit StringCopy( const char* s ) : ptr_(strdup(s)) { }
  ~StringCopy() { free(ptr_); }
};

bool 
ErrOptParser::parseOpts(const char* optarg)
{
  bool summaryOK = true;

  char** optRegister = new char*[optmap_.size()+1];
  int i=0;
  for(MapType::iterator j=optmap_.begin(); j != optmap_.end(); ++i,++j)
    optRegister[i] = const_cast<char*>(j->first.c_str());
  optRegister[i] = 0;    // getsubopt requires null at end


  // getsubopt() likes to modify the string it works on, so we have
  // to make a local copy -- which we own and must free.
  StringCopy optcopy(optarg);

  char* modopt = optcopy.ptr_; // copy of the pointer that getsubopt may modify
  char* value;

  while(*modopt){
    int c = getsubopt(&modopt, optRegister, &value);
    if(c == -1){
      cerr << "Unrecognized argument \"" << value 
		<< "\" for error option" << endl;
      summaryOK = false;
    }
    else {
      summaryOK = summaryOK && processOption(optRegister[c], value);
    }
  }

  delete[] optRegister;
  return summaryOK;
}

bool 
ErrOptParser::processOption( const std::string& name, const char* value)
{
  MapType::iterator opti = optmap_.find(name);
  
  bool stat;
  if(opti == optmap_.end()){
    cerr << "Internal error: unregistered suboption \""
	      << name << "\" encountered" << endl;
    stat = false;
    assert(opti != optmap_.end());
  }
  else{
    stat = opti->second->processValue(value);
    if(!stat)
      cerr << "Error converting \"" << value
		<< "\" to value for option \"" << name << "\"" << endl;
  }
  return stat;
}
  
template<class T> 
bool 
ErrOptParser::OptHelper<T>::processValue(const char* value)
{
  if(!value)
    switch(optref_.type()){
    case ErrOpt::Required:
      cerr << "Missing value for error suboption \"" 
	   << optref_.name() << "\"" << endl;
      return false;
    case ErrOpt::Optional:
    case ErrOpt::Prohibited:
      ErrOptParser::setPresent( optref_, true );
      return true;
    }
  else{
    istringstream vis(value);
    switch(optref_.type()){
    case ErrOpt::Required:
    case ErrOpt::Optional:
      {
	vis >> optref_.value();
	// The following test is required to get consistent behaviors on
	// all platforms.  We want to require not only that a value of the
	// given type could be constructed from the string, but also that
	// _all_ the string was consumed to do so.  Ideally we would require
	//
	//   vis.rdstate() == ios::eofbit
	// 
	// thus checking for no actual error but requiring EOF.
	// Unfortunately, for some reason, for the input of integers from 
	// stringstreams, the stream is often not left with eofbit set even when
	// all its characters have been consumed.  So we have to use the
	// following construct to check for 1) no error, and 2) EOF.
	bool converted = vis && vis.peek() == EOF;
	ErrOptParser::setPresent( optref_, converted );
	if ( !converted ) {
	  cerr << "Error parsing value \"" << value << "\" for suboption \""
	       << optref_.name() << "\"" << endl;
	}
	return converted;
      }
    case ErrOpt::Prohibited:
      cerr << "Error suboption \"" 
	   << optref_.name() << "\" does not accept a value" << endl;
      return false;
    }
  }
}

bool 
ErrOptParser::BoolHelper::processValue(const char* value)
{
  if (!value) {
    // There's no value associated with the suboption in the input string.
    if ( negated_ ) {
      // That's exactly right for the negated suboption.
      ErrOptParser::setPresent( optref_, true );
      optref_.value() = false;
      return true;
    }
    else {
      switch (optref_.type()) {
      case ErrOpt::Required:
	cerr << "Missing value for error suboption \"" 
	     << optref_.name() << "\"" << endl;
	return false;
      case ErrOpt::Optional:
      case ErrOpt::Prohibited:
	ErrOptParser::setPresent( optref_, true );
	optref_.value() = true;
	return true;
      }
    }
  }
  else {
    // There is a value.  
    if ( negated_ ) {
      // A value is not permitted for the negated "noFoo" form of a
      // negatable option, because double negatives (e.g. "noFoo=false")
      // would confuse users.
      cerr << "Negated Boolean suboption \"no" << optref_.name() 
	   << "\" does not accept a value" << endl;
      return false;
    }
    else {
      switch (optref_.type()) {
      case ErrOpt::Required:
      case ErrOpt::Optional:
	{
	  bool converted = getbool(value, optref_.value());
	  ErrOptParser::setPresent( optref_, converted );
	  if ( !converted ) {
	    cerr << "Error parsing value \"" << value << "\" for suboption \""
		 << optref_.name() << "\"" << endl;
	  }
	  return converted;
	}
      case ErrOpt::Prohibited:
	cerr << "Error suboption \"" << optref_.name() 
	     << "\" does not accept a value" << endl;
	return false;
      }
    }
  }

  assert(0);
  return false;
}

bool
ErrOptParser::BoolHelper::getbool( const std::string& s, bool& b)
{
  if ( s == "true" || s == "1" ) {
    b = true;
    return true;
  }
  else if ( s == "false" || s == "0" ) {
    b = false;
    return true;
  }

  return false;
}
