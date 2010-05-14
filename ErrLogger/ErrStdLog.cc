//--------------------------------------------------------------------------
//
// Description:
//
//  Error & logging standard implementation
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Scott Metzler
//
//      $Id: ErrStdLog.cc,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
//
// Copyright Information:
//	Copyright (C) 1998, 2002   California Institute of Technology
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "ErrLogger/ErrStdLog.hh"

//-----------------
// C/C++ Headers --
//-----------------
#include <iostream>
#include <sstream>
#include <assert.h>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "ErrLogger/ErrStreamWrapper.hh"
using std::cerr;
using std::cout;
using std::endl;
using std::ios;
using std::ostream;
using std::ostringstream;
using std::stringstream;


//		--------------------------------
// 		-- Local Variable Definitions --
//		--------------------------------
static ostringstream devnull;

//		---------------------------------
// 		-- Global Variable Definitions --
//		---------------------------------

//		---------------------------------
// 		-- Global Function Definitions --
//		---------------------------------

ErrStdLog::ErrStdLog( Severity theSeverity )
  : ErrLog( theSeverity )
  , _myStreamCout( new ErrStreamWrapper(cout,this) )
  , _myStreamCerr( new ErrStreamWrapper(cerr,this) )
  , _myDevnull( new ErrStreamWrapper(devnull,this) )
{
  // Establish this as the active implementation.
  ErrLog::_setImplementation( this );
}


ErrStdLog::ErrStdLog( ostream& theStream, Severity theSeverity )
  : ErrLog( theSeverity )
  , _myStreamCout( new ErrStreamWrapper(theStream,this) )
  , _myStreamCerr( new ErrStreamWrapper(theStream,this) )
  , _myDevnull( new ErrStreamWrapper(devnull,this) )
{
  // Establish this as the active implementation.
  ErrLog::_setImplementation( this );
}


ErrStdLog::~ErrStdLog()
{
  // There can be only one, and it should be this one.
  assert( ErrLog::_implementation() == this );

  // Remove this one, leaving the default implementation in place
  // (in case subsequently called destructors issue ErrMsgs).
  ErrLog::_setImplementation(0);

  // cleanup this class
  delete _myStreamCout;
  delete _myStreamCerr;
  delete _myDevnull;
}

ErrStream& 
ErrStdLog::doMsg( ErrLog::Severity severity, 
		  const char* facility, 
		  int code )
{
  // initialize result
  ErrStream* resultPtr = 0;

  // select which stream to write to
  if ( logging(severity,facility,code) ) {
    // pick the correct stream
    if (severity>=warning)
      resultPtr = _myStreamCerr;
    else
      resultPtr = _myStreamCout;
  } // end if logging
  else { // not logging
    // send to local char stream after resetting to start over
    devnull.rdbuf()->pubseekoff(0, ios::beg,ios::out);
    resultPtr = _myDevnull;
  } // end if

  // check that a result was set
  assert( resultPtr!=0 );

  // return result
  return *resultPtr;
}

void 
ErrStdLog::doEndmsg( const char* text, ErrStream& stream )
{
  // Since the streams defined by this class are all ErrStreamWrapper,
  // they do output as-you-go and should never pass text back to here.
  assert( text == 0 );

  // But we can handle it if assertions are compiled away...
  if ( text != 0 ) {
    static_cast<ostream&>(stream) << text;
  }

  // Terminate all messages with endl (at least for now; eventually we
  // might decide to -- optionally? -- terminate low-severity messages
  // with just '\n' for improved performance).
  static_cast<ostream&>(stream) << endl;

  // Apply the default ErrLog termination policy.
  ErrLog::_defaultTermination( stream.severity(), stream.code() );
}
