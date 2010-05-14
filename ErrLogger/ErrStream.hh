#ifndef ERRSTREAM_HH
#define ERRSTREAM_HH

//------------------------------------------------------------------------
//
// Description:
//	Class ErrStream :
//	This class encapsulates the ostream class and provives a
//	way to "trap" messages by using "endmsg".  This class is
//	sub-classes by OstErrLogger/OstLogStream.  When an "endmsg"
//	is processed, the message is considered finished.  "endmsg"
//	is distinct from "endl" and '\n'. 
//
// **********************************************************************
//
// Proper use:
//	ErrStream should not be used directly.  It should be
//	accessed through the ErrLog interface.
//
//	ErrMsg(warning) << "some resource not available" << endmsg;
//
//	It is important to end all error messages with "endmsg".
//	If a user does not finish an error message with "endmsg",
//	a warning will be printed that looks like:
//
//	 *** ErrLogger Warning *** 
//	The ErrStream was previously used without
//	endmsg being called to end the error message
//	Please see ErrLogger/ErrStream.hh for more
//	information.
//
//	 Error Location: ErrLogTest.cc::47
//	 *** ErrLogger Warning Finished ***
//
//	If you see the above warning, you need to end your error
//	message with << endmsg;, as described above.
//
// ***********************************************************************
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Scott Metzler		Original author
//      $Id: ErrStream.hh,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
//
// Copyright Information:
//	Copyright (C) 1998-2004	California Institute of Technology
//
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------
#include <string>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "ErrLogger/ErrLog.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include <iosfwd>


//		---------------------
// 		-- Class Interface --
//		---------------------

class ErrStream {
protected:

  ErrStream( std::ostream*, ErrLog* );
  ErrStream( const ErrStream& );

public:
  virtual ~ErrStream() {}

  void activate() {_self = this;}
    // set this ErrStream object as the current one

  // conversion operator
  operator std::ostream&() { return *_myOstream; }

  ErrLog::Severity severity() const { return _severity; }
  const char*      facility() const { return _facility.c_str(); }
  int              code()     const { return _code; }

protected:

  static ErrStream* activeStream() { return _self; }
  
  // virtual methods
  virtual void doEndmsg( ) = 0;

  // helper functions

  virtual void setStream( std::ostream* os ) { _myOstream = os; }

  bool inUse() const { return _inUse; }
  void setInUse(bool use) { _inUse = use; }

  ErrLog* logger() { return _myLogger; }
  void loggerEndmsg( const char* text ) { logger()->doEndmsg( text, *this ); }

  // Set message parameters to keep around for endmsg.  The facility string 
  // is copied to avoid pointer validity issues.
  void setParameters( ErrLog::Severity severity, 
		      const char* facility, 
		      int code );
  
  // This is a pass-through function for subclasses to use so that they
  // can benefit from ErrStream's friendship with ErrLog.
  void defaultTermination() 
  { ErrLog::_defaultTermination( severity(), code() ); }

  // friend relationships needed to implement "endmsg"
  friend std::ostream& operator<<( std::ostream& os, void (*fp)(ErrStream& es) );
  friend void endmsg( ErrStream& es );
  friend class ErrLog;

  // static variables
  static ErrStream* _self;
  static bool _inUse;

  // member variables
  std::ostream* _myOstream;
  ErrLog*  _myLogger;

  ErrLog::Severity _severity;
  std::string      _facility;
  int              _code;

private:

  // Disable assignment operator.
  ErrStream& operator=( const ErrStream& );

};

#endif // ERRSTREAM_HH
