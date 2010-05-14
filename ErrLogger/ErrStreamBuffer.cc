//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: ErrStreamBuffer.cc,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
//
// Description:
//	ErrStreamBuffer implementation
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
#include "ErrLogger/ErrStreamBuffer.hh"

//-----------------
// C/C++ Headers --
//-----------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "ErrLogger/ErrLog.hh"
using std::ends;

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
ErrStreamBuffer::ErrStreamBuffer( ErrLog* logger )
  : ErrStream( 0, logger )
{
  // Initialization-order constraints prohibited establishing the
  // connection between the base class stream pointer and the ostringstream
  // above.  Instead we set it now that we know the ostringstream is initialized.
  setStream( &_bufferStream );
}

//--------------
// Destructor --
//--------------
ErrStreamBuffer::~ErrStreamBuffer()
{ }

//-------------
// Methods   --
//-------------
    
void
ErrStreamBuffer::doEndmsg()
{
  if ( logger() != 0 ) {
    // We require that we be tied back to a real ErrLog implementation.
    // (If not, we do nothing except reset the buffer.)

    // Make sure the buffer is null-terminated.
    _bufferStream << ends;

    // Extract the string from the buffer and hand it to the ErrLog.
    // Note that this freezes the buffer and relinquishes ownership of the
    // string.
    const char* text = _bufferStream.str().c_str();
    loggerEndmsg( text );

    // Unfreeze the buffer.
    // _bufferStream.rdbuf()->freeze(false) ;
  }

  // Reset the buffer for the next message.
  _bufferStream.clear() ;
  _bufferStream.seekp(0) ;

}
