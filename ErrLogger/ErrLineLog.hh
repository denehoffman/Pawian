#ifndef ERRLINELOG_HH
#define ERRLINELOG_HH
//--------------------------------------------------------------------------
// $Id: ErrLineLog.hh,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
//
// Description:
//	Class ErrLineLog :
//	Implementation of ErrLog based on ErrStdLog.  In addition to
//	allowing a user to set the severity of messages sent to the
//	output streams (cout and cerr), this version writes the facility
//      (file name and line number, if ErrMsg is used) and code (if non-zero)
//      in the same way that FrameLogger/AppErrLog does.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Michael Kelsey <kelsey@slac.stanford.edu>
//
// Copyright Information:
//	Copyright (C) 1999  Princeton University
//      Copyright (C) 2004  California Institute of Technology
//
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------

//----------------------
// Base Class Headers --
//----------------------
#include "ErrLogger/ErrStdLog.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class ErrStream; // declare


//		---------------------
// 		-- Class Interface --
//		---------------------

class ErrLineLog : public ErrStdLog {
  public:
    // constructor
    ErrLineLog( ErrLog::Severity theSeverity = ErrLog::warning );

    // destructor
    virtual ~ErrLineLog();

  protected:

    virtual ErrStream& doMsg( ErrLog::Severity severity, 
			      const char* facility, 
			      int code );

  // All other interface uses ErrStdLog implementation

  private:
    // dis-allow copy construction and assignment
    ErrLineLog( const ErrLineLog& );
    ErrLineLog& operator=( const ErrLineLog& );
};

#endif /* ERRLINELOG_HH */
