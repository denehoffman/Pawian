//--------------------------------------------------------------------------
// $Id: ErrLineLog.cc,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
//
// Description:
//
//  Error & logging implementation which writes facility/code to messages.
//  Based on ErrLogger/ErrStdLog.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Michael Kelsey <kelsey@slac.stanford.edu>
//
// Copyright Information:
//	Copyright (C) 1999  Princeton University
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "ErrLogger/ErrLineLog.hh"

//-----------------
// C/C++ Headers --
//-----------------
#include <assert.h>
#include <string.h>
#include <iostream>

//-------------------------------
// Collaborating Class Headers --
//------------------------------
#include "ErrLogger/ErrStream.hh"
using std::ostream;


// --------------------------------
// -- Constructor and Destructor --
// --------------------------------
ErrLineLog::ErrLineLog( Severity theSeverity )
  : ErrStdLog(theSeverity) 
{ }		// No actions beyond base class


ErrLineLog::~ErrLineLog() 
{
  // There can be only one, and it should be this one.
  assert( ErrLog::_implementation() == this );
}

//-------------
//-- Actions --
//-------------
ErrStream& 
ErrLineLog::doMsg( ErrLog::Severity severity, 
		   const char* facility, 
		   int code )
{
  // Use base class to select stream, then write info
  ErrStream& theStream = ErrStdLog::doMsg(severity,facility,code);

  ostream& toWrite( theStream );	// Cast to ostream for writing;

  // facility is the filename with the long path removed
  char shortName[100];
  char *fromLastSlash = (char*) strrchr(facility, '/');
  if(fromLastSlash != NULL)
    {
      strncpy(shortName, &fromLastSlash[1], 100);
    }
  else
    {
      strncpy(shortName, facility, 100);
    }

  toWrite << shortName;
  if ( code != 0 ) {
    toWrite << " [" << code << "]";
  }
  toWrite << ":";

  return theStream;
}
