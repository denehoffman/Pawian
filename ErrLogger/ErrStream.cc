//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: ErrStream.cc,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
//
// Description:
//	Template implementation
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Scott Metzler (originator)
//	Gregory Dubois-Felsmann  (extensive revisions 2004)
//
// Copyright Information:
//	Copyright (C) 1998-2004   California Institute of Technology
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "ErrLogger/ErrStream.hh"

//-----------------
// C/C++ Headers --
//-----------------
#include <iostream>
#include <stdlib.h>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "ErrLogger/ErrLog.hh"
using std::ostream;

//		-----------------------------------------------
// 		-- Static Data & Function Member Definitions --
//		-----------------------------------------------

ErrStream* ErrStream::_self = 0;
bool       ErrStream::_inUse = false;

//		----------------------------------------------------
// 		-- Public & Protected Function Member Definitions --
//		----------------------------------------------------

ErrStream::ErrStream( ostream* os, ErrLog* logger )
  : _myOstream(os)
  , _myLogger(logger)
  , _severity(ErrLog::debugging)
  , _facility()
  , _code(0)
{
}

ErrStream::ErrStream( const ErrStream& rhs )
  : _myOstream(rhs._myOstream)
  , _myLogger(rhs._myLogger)
  , _severity(rhs._severity)
  , _facility(rhs._facility)
  , _code(rhs._code)
{
}

void
ErrStream::setParameters( ErrLog::Severity severity,
			  const char*      facility,
			  int              code )
{
  _severity = severity;
  _facility = facility;
  _code     = code;
}
 

// Implementation of the global operator<< that allows "endmsg" to be used
// with a plain ostream as a manipulator.  It relies on there being only one
// active ErrStream at any given time.

ostream& 
operator<<( ostream& os, void (* fp)(ErrStream&) )
{
  (*fp)(*ErrStream::_self);
  return *ErrStream::_self; 
}

// Implementation of the global function that is actually called when the
// endmsg manipulator is used.

void endmsg( ErrStream& es )
{
  es.doEndmsg();
  es.setInUse(false);
}
