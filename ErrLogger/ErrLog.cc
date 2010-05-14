//--------------------------------------------------------------------------
//
// Description:
//
//  Error & logging abstract base class
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Bob Jacobsen
//      $Id: ErrLog.cc,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "ErrLogger/ErrLog.hh"

//-----------------
// C/C++ Headers --
//-----------------
#include <string.h>
#include <assert.h>

#include <iostream>
#include <sstream>
#include <string>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "ErrLogger/ErrStreamWrapper.hh"
#include "ErrLogger/ErrOpt.hh"
#include "ErrLogger/ErrOptParser.hh"
using std::cerr;
using std::cout;
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


// static size_t
// strnlen( const char* string, size_t maximum )
// {
//   for ( size_t length = 0; length < maximum; length++ ) {
//     if (string[length] == '\0') return length;
//   }
//   return maximum;
// }



ErrLog::ErrLog( ErrLog::Severity minSeverity )
  : _useStdio(true)
  , _useStdioOpt( new ErrOptRef<bool>("stdio", ErrOpt::Optional, _useStdio) )
  , _itsErrOptParser( new ErrOptParser )
{
  // It is assumed that any subclass constructor calling this protected
  // constructor will be establishing itself as the active ErrLog instance.
  // So we allow it to set the default minimum severity in the static front
  // end to ErrLog.
  _minSeverity = minSeverity;
}


ErrLog::~ErrLog()
{
  delete _itsErrOptParser;
  delete _useStdioOpt;
}


ostream& 
ErrLog::msg( ErrLog::Severity severity, 
	     const char* facility, 
	     int code )
{
  if ( ErrStream::activeStream() != 0 
       && ErrStream::activeStream()->inUse() ) {
    // A new message is being started before the last one has been
    // completed with endmsg.

    // 1) Flush what has come of the original message so far.
    // 2) Optionally write a warning about this.
    // 3) Continue on to establish the new message.

    // 1) We complete the previous pending message.  This is important so
    //    that, if it's buffered, it will be written with its original
    //    severity and other parameters.
    ErrStream* prevStream = ErrStream::activeStream();

    if ( _nestedWarnings ) {
      static_cast<ostream&>(*prevStream) 
	<< " [interrupted by nested message from " << facility << "]";
    }
    endmsg( *prevStream );

    if ( _nestedWarnings ) {
      // 2) Warn about this usage.
      
      // Need to trap the facility and code, because they will be 
      // overwritten by the ErrMsg call itself.
      std::string prevFacility( prevStream->facility() );
      int prevCode = prevStream->code();

      ErrMsg(warning) << "\n *** ErrLogger Warning *** \n"
		      << "The ErrStream was previously used without\n"
		      << "endmsg being called to end the error message.\n"
		      << "Please see ErrLogger/ErrStream.hh for more\n"
		      << "information.\n\n"
		      << " Error Location/Facility: " << prevFacility
		      << " [code " << prevCode << "]\n"
		      << " *** ErrLogger Warning Finished ***\n"
		      << endmsg;
    }

    // 3) Now we can proceed.
  }

  // initialize result
  ErrStream* resultPtr = 0;

  if (_implementation() != 0) {
    resultPtr = & ( _implementation()->doMsg( severity, facility, code ) );
  }
  else {
    resultPtr = _getDefaultStream( severity, facility, code );
  } // end else

  // check that resultPtr is set
  assert( resultPtr!=0 );

  // Save the message parameters.
  _setParameters( *resultPtr, severity, facility, code );

  resultPtr->activate();

  // mark this stream as in use
  resultPtr->setInUse(true);
      
  return static_cast<ostream&>( *resultPtr );
}



ostream& 
ErrLog::msg( ErrLog::Severity severity, 
	     const char* facility, 
	     const char* code )
{
  // This implementation is provided for backward compatibility with the
  // old, misunderstanding-derived notion that the "const char* code"
  // form of msg() should be the primary one.  Now (June 2002) that we
  // have switched to the original intent, where "int code" is the 
  // primary interface, we still provide this interface, in a way which
  // does almost exactly what the old one did in practice in all the
  // existing ErrLog implementations: append "(" + code + ")" to the
  // facility.

  std::string newfac( facility );
  newfac += "(";
  newfac += code;
  newfac += ")";

  return ErrLog::msg( severity, newfac.c_str(), 0 );
}


bool
ErrLog::logging( ErrLog::Severity severity, 
		 const char* facility, 
		 const char* code )
{
  // See msg( ErrLog::Severity severity, const char* facility, 
  //          const char* code ) above.  Note that these string
  // manipulations are a lot of work to do if you just want to check 
  // to see whether a particular class of message is being logged!

  std::string newfac( facility );
  newfac += "(";
  newfac += code;
  newfac += ")";

  return ErrLog::logging( severity, newfac.c_str(), 0 );
}



const char*
ErrLog::severityName( ErrLog::Severity severity )
{
  switch ( severity ) {
  case ErrLog::debugging:
    return "debugging";
  case ErrLog::trace:
    return "trace";
  case ErrLog::routine:
    return "routine";
  case ErrLog::warning:
    return "warning";
  case ErrLog::error:
    return "error";
  case ErrLog::severe:
    return "severe";
  case ErrLog::fatal:
    return "fatal";
  }

  return "?invalid-severity?";
}



char
ErrLog::severityLetter( ErrLog::Severity severity )
{
  switch ( severity ) {
  case ErrLog::debugging:
    return 'D';
  case ErrLog::trace:
    return 'T';
  case ErrLog::routine:
    return 'R';
  case ErrLog::warning:
    return 'W';
  case ErrLog::error:
    return 'E';
  case ErrLog::severe:
    return 'S';
  case ErrLog::fatal:
    return 'F';
  }

  return '?';
}



bool 
ErrLog::doLogging( ErrLog::Severity severity,
		   const char* facility,
		   int code )
{
  return ErrLog::_getDefaultLogging( severity, facility, code );
}


ErrStream* 
ErrLog::_getDefaultStream( ErrLog::Severity severity,
			   const char* facility,
			   int code )
{
  // initialize result
  ErrStream* result = 0;

  // here use the default implementation
  if ( _myStreamCerr==0 )
    _myStreamCerr = new ErrStreamWrapper(cerr);
  if ( _myStreamCout==0 )
    _myStreamCout = new ErrStreamWrapper(cout);
  if ( _myDevnull==0 )
    _myDevnull    = new ErrStreamWrapper(devnull);

  // select which stream to write to
  if ( logging( severity, facility, code ) ) {
    // pick the correct stream
    if ( severity>=warning )
      result = _myStreamCerr;
    else
      result = _myStreamCout;
  } // end if logging
  else { // not logging
    // send to local char stream after resetting to start over
    devnull.rdbuf()->pubseekoff(0, ios::beg,ios::out);
    result = _myDevnull;
  } // end if

  return result;
}


void
ErrLog::_setParameters( ErrStream& es,
			ErrLog::Severity severity,
			const char* facility,
			int code )
{
  es.setParameters( severity, facility, code );
}

void
ErrLog::_defaultTermination( ErrLog::Severity severity, 
			     int code )
{
  _executeDefaultTermination( _getDefaultTermination( severity, code ) );
}

int
ErrLog::_getDefaultTermination( ErrLog::Severity severity, 
				int code )
{
  int result = 0;
  if ( severity == ErrLog::severe ) {
    // Default behavior for this level is to exit with a non-zero status.
    // That status is determined based on the "code".  For a code of zero,
    // the default exit code is used.
    if ( code == 0 ) {
      result = DefaultExitCode;
    }
    else {
      // Otherwise we take the low-order seven bits of the "code", so
      // that we stay in the range shells understand as exit()-driven
      // process terminations.
      result = code & 0x7f;

      // Except we never allow that to be zero, either.  We assign a
      // special code as a warning.
      if ( result == 0 ) result = DefaultExitCode - 1;
    }
  }
  else if ( severity == ErrLog::fatal ) {
    // Return the value that indicates "abort".
    result = -1;
  }

  return result;
}

void
ErrLog::_executeDefaultTermination( int value )
{
  if ( value == -1 ) ::abort();
  else if ( value > 0 ) ::exit( value );
}

void
ErrLog::turnOffWarnings()
{
  _nestedWarnings = false;
}

bool
ErrLog::_configureOptions()
{
  bool ok = true;
  ok = ok && _errOptParser().registerNegatableOption( *_useStdioOpt );
  return ok;
}

bool 
ErrLog::parseOpts( const char* optarg ) 
{
  return _errOptParser().parseOpts(optarg);
}


// Protected static member initializers.

ErrStream*       ErrLog::_myStreamCerr = 0;
ErrStream*       ErrLog::_myStreamCout = 0;
ErrStream*       ErrLog::_myDevnull = 0;
ErrLog::Severity ErrLog::_minSeverity = ErrLog::debugging;


// Private static member initializers.

ErrLog*          ErrLog::_implInstance   = 0;
bool             ErrLog::_nestedWarnings = true;

