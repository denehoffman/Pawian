#ifndef ERRSTREAMBUFFER_HH
#define ERRSTREAMBUFFER_HH

//--------------------------------------------------------------------------
//
// $Id: ErrStreamBuffer.hh,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
//
// Environment:
//      This software was developed for the BaBar collaboration.  If you
//      use all or part of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 2004   California Institute of Technology
//
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------
#include <sstream>

//----------------------
// Base Class Headers --
//----------------------
#include "ErrLogger/ErrStream.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class ErrLog;


//		---------------------
// 		-- Class Interface --
//		---------------------

/**
 *  This class implements a stream that buffers output when written to,
 *  making it available back to the ErrLog implementation when a message is
 *  completed with endmsg.
 *
 *  It is used in the special message logging schemes in OstErrLogger and
 *  QtrErrLog.
 */

class ErrStreamBuffer : public ErrStream {
public:

  ErrStreamBuffer( ErrLog* );
  virtual ~ErrStreamBuffer();

protected:

  // virtual methods
  virtual void doEndmsg( );

private:

  std::ostringstream _bufferStream;

  // Disable copy-construction and assignment.
  ErrStreamBuffer( const ErrStreamBuffer& );
  ErrStreamBuffer& operator=( const ErrStreamBuffer& );

};

#endif // ERRSTREAMBUFFER_HH
