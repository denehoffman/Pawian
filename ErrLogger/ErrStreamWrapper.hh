#ifndef ERRSTREAMWRAPPER_HH
#define ERRSTREAMWRAPPER_HH

//--------------------------------------------------------------------------
//
// $Id: ErrStreamWrapper.hh,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
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
#include <iosfwd>
class ErrLog;


//		---------------------
// 		-- Class Interface --
//		---------------------

/**
 *  This class implements a wrapper for a stream that actually produces
 *  output when written to -- an std::ofstream or std::cout/std::cerr/clog.
 *
 *  It is used in the default ErrLog front-end as well as in the std::cout/cerr
 *  flavors of ErrLog.
 */

class ErrStreamWrapper : public ErrStream {
public:

  ErrStreamWrapper( std::ostream&, ErrLog* = 0 );
  virtual ~ErrStreamWrapper();

protected:

  // virtual methods
  virtual void doEndmsg( );

private:

  // Disable copy-construction and assignment.
  ErrStreamWrapper( const ErrStreamWrapper& );
  ErrStreamWrapper& operator=( const ErrStreamWrapper& );

};

#endif // ERRSTREAMWRAPPER_HH
