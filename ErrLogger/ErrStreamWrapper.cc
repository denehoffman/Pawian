//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: ErrStreamWrapper.cc,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
//
// Description:
//	ErrStreamWrapper implementation
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Gregory Dubois-Felsmann  (originator)
//
// Copyright Information:
//	Copyright (C) 2004   California Institute of Technology
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "ErrLogger/ErrStreamWrapper.hh"

//-----------------
// C/C++ Headers --
//-----------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "ErrLogger/ErrLog.hh"
using std::endl;
using std::ostream;

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
ErrStreamWrapper::ErrStreamWrapper( ostream& os, ErrLog* logger )
  : ErrStream( &os, logger )
{
}

//--------------
// Destructor --
//--------------
ErrStreamWrapper::~ErrStreamWrapper()
{ }

//-------------
// Methods   --
//-------------
    
void
ErrStreamWrapper::doEndmsg()
{
  if ( logger() != 0 ) {
    // We are tied back to a real ErrLog implementation.
    // Because this is a wrapper for a real-output stream, there is no
    // additional text to output: hence the null pointer passed through.
    // The ErrLog is responsible for adding a newline at the end of the
    // message, and so can choose between endl and '\n'.
    loggerEndmsg( 0 );
  } 
  else {
    // No ErrLog is available to do the wrap-up.

    // Write a newline, using endl to be conservative about buffering.
    static_cast<ostream&>(*this) << endl;

    // Handle the termination.  This has to be done here, again because there
    // is no stateful ErrLog instance to go back to for the termination
    // decision for the current message.
    defaultTermination();
  }
}
