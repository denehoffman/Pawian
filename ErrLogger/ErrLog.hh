#ifndef ERRLOG_HH
#define ERRLOG_HH

//--------------------------------------------------------------------------
// Description:
//	Class ErrLogger:
//		Interface description for the abstract logging interface
//
// BaBar Usage:
//
//	BaBar policy for using ErrLogger is to interact with it only
//	through the ErrMsg macro.
//
//	ErrMsg(severity) << "the message" << endmsg;
//
//	The error message must be terminated by endmsg as above.
//
//	The ErrMsg macro returns an ostream& so that it is possible
//	to include any constructs that can be sent to cout or cerr.
//	For example, the string classes and user-defined operator<<
//	are available to users while using ErrMsg.
//
// General Usage:
//
//	Available constructs are a logging call:
//
//		ErrLog::msg(severity, facility, code) << "Foo" << endmsg;
//
//	and a test to see if a particular logging call will do anything:
//
//		if (ErrLog::logging(severity, facility, code)) ...
//
//	All error messages should be completed with "endmsg".  For multiple-
//	line messages use '\n' for line breaks in preference to endl.  The 
//      endmsg command will supply a final "endl" for error loggers for which
//      it's appropriate (e.g., if writing to the console).  Users should in 
//      particular NOT terminate messages with "some text" << endl << endmsg;.
//
//
//	"severity" is a member of an enum; see the definition below
//	"facility" names the source of the message; in the ErrMsg() form,
//                 this is set to be the filename and line number
//	"code"     an integer distinguishing the individual errors possible;
//	           no meanings are predefined at this level
//
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Bob Jacobsen	(original author)
//	Scott Metzler	(addition of ErrStream)
//	Gregory Dubois-Felsmann   (ErrOpt; conversion of "code" to int)
//      $Id: ErrLog.hh,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
//
// Copyright Information:
//	Copyright (C) 1998  Lawrence Berkeley National Laboratory
//	Copyright (C) 1998-2004  California Institute of Technology
//
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------
#include <assert.h>
#include <iostream>

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class ErrStream;
class ErrOptParser;
template <class T> class ErrOptRef;


// Use the ErrMsg and ErrLogging macros for interacting with the error logger
// in a syntactically economical way.

// The next two macros allow the substitution of the value of the 
// __LINE__ macro as a quoted character string instead of an integer.
// As odd as this may seem, it's a well-established ISO C trick.

// 1) replaced by its argument, treated as a quoted string, without any
//    nested expansions applied (they are suppressed by the "#" operator).
#define ERRLINE_HACK_1(line)   #line

// 2) allows its argument (which will be "__LINE__") to have macro expansion
//    performed on it.  Hence HACK_2(__LINE__) becomes, e.g., HACK_1(76).
#define ERRLINE_HACK_2(line)   ERRLINE_HACK_1(line)


#ifdef ErrMsg
#undef ErrMsg
#endif
#define ErrMsg(sev) ErrLog::msg( ErrLog::sev, __FILE__ "(" ERRLINE_HACK_2(__LINE__) ")", 0 )

#ifdef ErrLogging
#undef ErrLogging
#endif
#define ErrLogging(sev) ErrLog::logging( ErrLog::sev, __FILE__ "(" ERRLINE_HACK_2(__LINE__) ")", 0 )

// [June 2002 note]  The old versions of these macros expanded to 
// ( ErrLog::sev, __FILE__, __LINE__ ) and so invoked the (..., int) forms
// of msg() and logging().  But at run time these forms were immediately
// vectored, by converting the code to a string, into the (..., const char*)
// forms, in which then in practice in all relevant ErrLog implementations the
// "code" string, containing the line number, was appended to the facility.
//
// The new forms of these macros preserve this behavior -- the line number
// is concatenated to the file name -- on the theory that the line number
// makes a silly code number and that the file+line combination is 
// reasonable to treat as a unit.


// The actual class

class ErrLog {
public:

  /**
   *  Message severity levels.
   */
  enum Severity { debugging=-1, 
		  trace=0, 
		  routine, 
		  warning, 
		  error,
		  severe,
		  fatal, 
		  MinSeverity=debugging,
		  MaxSeverity=fatal };

  //	fatal:		A condition exists that casts doubt on the integrity
  //                    of processing and the reliability of any results.
  //                    No recovery strategy can be attempted, because the
  //                    state of the program may be indeterminate.
  //
  //                    It is intended that "fatal" errors be reported
  //                    primarily when it is likely that the condition was
  //                    caused by incorrect code, or when the root cause cannot
  //                    be determined.  Dealing with the error is likely to
  //                    require debugging or other expert investigation.
  // 
  //                    In most implementations, this will be configured to
  //                    result in application termination with ::abort() 
  //                    immediately after the error message is completed with
  //                    "endmsg", producing a core dump usable to obtain a
  //                    traceback.
  //
  //			Programmers should generally not call ::abort() or
  //                    ::exit() themselves following the message.

  //    severe:         A condition exists such that the general ability of
  //                    the program to continue work is compromised: the
  //                    present computation and most or all subsequent ones
  //                    are likely to fail or produce erroneous results.
  //
  //                    It is intended that "severe" errors be reported 
  //                    primarily when external conditions (command line
  //                    arguments, control files, data files, the contents of
  //                    databases) are identified to be incompatible with
  //                    proceeding.  The programmer should ensure that the
  //                    message written is as informative as possible.
  //
  //                    In most implementations, this will be configured to
  //                    result in application termination immediately after
  //                    the error message is completed with "endmsg" and
  //                    delivered, with a non-zero exit status.  The ErrLog
  //                    implementation is free to determine the exit status
  //                    to use, from the range 1-126.  It may take the "code"
  //                    value from the ErrLog::msg( , , code ) into account.
  //
  //			Programmers should generally not call ::abort() or
  //                    ::exit() themselves following the message.

  //	error:		A condition exists such that requested result or
  //			action can not be produced.  However, future 
  //                    processing may be able to succeed.

  //	warning:	The result is produced, but may not be
  //			what's desired due to an unexpected condition

  //	routine:	Nothing known to be wrong with the result;
  //			messages that are always produced in normal
  //			operation

  //	trace:		Messages about the flow of program control
  //			and which optional operations took place.
  //			(This is the default if nothing is defined)
  
  //	debugging:	Information in addition to the above


  /**
   *  This constant is the exit code to be used when a message severity
   *  is such that an exit(!=0) -- as opposed to an abort() -- is merited,
   *  but no message-specific exit code can be determined from the 
   *  message parameters.
   */
  enum { DefaultExitCode=125 };

//-----------------------------------------------------------------------
// Static public member functions ---------------------------------------
//-----------------------------------------------------------------------

public:

  // These member functions are the public interface to error logging.

  /**
   *  Returns a stream to which to write a message, with the parameters
   *  given.  The message must be terminated later by performing the
   *  "<< endmsg" operation on the resulting stream.  Actual output may
   *  not occur until the endmsg is executed.  The meaning of the "facility"
   *  and "code" parameters is not defined by the ErrLog base class, nor is
   *  is guaranteed that their values will be reflected in any output
   *  generated for the message.
   *
   *  The "facility" pointer must be guaranteed by the caller to remain
   *  valid until the message is completed with "endmsg".  (I.e., this
   *  function is not guaranteed to make a copy of the string.)
   */
  static std::ostream& msg( ErrLog::Severity severity, 
		       const char* facility, 
		       int code );

  /** 
   *  DEPRECATED
   *  Backward-compatibility function for an old API.
   *  Calls through to msg( severity, facilty + code, 0 )
   *  Do not use in new code.
   */
  static std::ostream& msg( ErrLog::Severity severity, 
			    const char* facility, 
			    const char* code );

  /**
   *  Tests whether actual logging output is enabled for a specific set
   *  of message parameters.  This function should be used before sending
   *  complex output to the result of ErrMsg() or ErrLog::msg(), to avoid
   *  the cost of executing the stream insertion operations (binary-to-
   *  string conversions) for the message, only to discard their output.
   */
  inline
  static bool logging( ErrLog::Severity severity, 
		       const char* facility, 
		       int code );

  /** 
   *  DEPRECATED
   *  Backward-compatibility function for an old API.
   *  Calls through to logging( severity, facilty + code, 0 )
   *  Do not use in new code.
   */
  static bool logging( ErrLog::Severity severity, 
		       const char* facility, 
		       const char* code );


  // The following are options controlling the global configuration of
  // the static "front end" to error logging.

  /**
   *  Disables the automatic warning message generated when a new ErrLog::msg()
   *  call is made before the previous one has been terminated with "endmsg".
   *  The use of this option is *strongly* discouraged; the correct response
   *  to such warnings is to fix the nested messages.
   */
  static void turnOffWarnings();


  // The following are utility functions that are probably of interest only
  // implementations (i.e., could be protected) but expose no internals and
  // are harmless, and maybe useful, to make public.

  /**
   *  Returns a static string corresponding to the enum member name of a
   *  severity level.  I.e., severityName(ErrLog::error) == "error".
   */
  static const char* severityName( ErrLog::Severity );

  /**
   *  Returns a single uppercase letter corresponding to the enum member name
   *  of a severity level.  I.e., severityName(ErrLog::error) == 'E'.
   */
  static char        severityLetter( ErrLog::Severity );



//-----------------------------------------------------------------------
// Instance public member functions -------------------------------------
//-----------------------------------------------------------------------

public:

  /**
   *  Provides parsing of a comma-delimited list of error logger
   *  behavior control options, in getsubopt() format.  The effects of
   *  the parsing depend on the specific ErrLog implementation and are
   *  established by it at constructor time using the _errOptParser()
   *  and _configureOptions() methods.
   */ 
  bool parseOpts( const char* optarg );


//-----------------------------------------------------------------------
// Friendship declarations ----------------------------------------------
//-----------------------------------------------------------------------

  friend class ErrStream;


//-----------------------------------------------------------------------
// Instance protected member functions ----------------------------------
//-----------------------------------------------------------------------

protected:

  // These members provide a basic set of functions for dealing with a
  // specific instance of ErrLog or a subclass.  They are really intended
  // only for use in configuring the active instance just after it's been
  // created.

  /**
   *  Protected constructor, for use in initializing subclasses.
   *
   *  @param minSeverity Minimum severity of messages that will be logged
   *                     by the default implementation.
   */
  ErrLog( ErrLog::Severity minSeverity = ErrLog::debugging );

  /**
   *  A virtual destructor is provided, but note that the quasi-singleton
   *  design of this class cannot guarantee that the active ErrLog instance
   *  will be deleted at application termination.
   */
  virtual ~ErrLog();

  // Functions used with ErrOptParser.

  virtual bool _configureOptions();
  ErrOptParser& _errOptParser() { return *_itsErrOptParser; }


//-----------------------------------------------------------------------
// Instance protected virtual member functions --------------------------
//-----------------------------------------------------------------------

protected:

  // These members provide the polymorphic interface for customizing the
  // behavior of an ErrLogger implementation.

  /** 
   *  Establishes an std::ostream that is to be used for logging message text,
   *  to be terminated with "endmsg".  Pure virtual.
   */
  virtual ErrStream& doMsg( ErrLog::Severity severity, 
			    const char* facility, 
			    int code ) = 0;

  /**
   *  Callback from ErrStream when a message has been completed, allowing
   *  the logger to decide what to do.
   *
   *  @param text   The full formatted message text, if the ErrStream is a
   *                "buffering" stream, or 0 (null pointer) if it is a
   *                "wrapper" stream, with output actually appearing as 
   *                each stream-output (<<) operation happens.
   *  @param stream The ErrStream to which the message had been directed.
   *                It's provided so that end-of-message output (std::endl or
   *                something more) can be written if appropriate.
   */
  virtual void doEndmsg( const char* text, ErrStream& stream ) = 0;

  /**
   *  Tests whether a message will actually be logged for a particular
   *  set of message parameters.  A default implementation is provided,
   *  based only on a test of the severity w.r.t. a minimum level, a
   *  static member of the ErrLog base class.
   *
   *  The sort of thing a subclass might do here is apply a "verbosity"
   *  control in addition to the severity test, perhaps suppressing
   *  repeated messages with the same parameters.
   */
  virtual bool doLogging( ErrLog::Severity severity, 
			  const char* facility, 
			  int code );


//-----------------------------------------------------------------------
// Static protected member functions ------------------------------------
//-----------------------------------------------------------------------

protected:

  // These members provide tools for ErrLog implementations to use.

  /** 
   *  ErrLog is a form of singleton.  Polymorphism of its static methods is
   *  arranged by calling implementation instances through this pointer.
   */
  inline static ErrLog* _implementation() { return _implInstance; }

  /** 
   *  ErrLog is a form of singleton.  The constructor of any subclass must
   *  register itself as the active instance by calling this function.
   */
  inline static void _setImplementation( ErrLog* implementation );

  /**
   *  Get the appropriate default stream for a given severity level.
   *  The facility and code arguments are not used.
   */
  static ErrStream* _getDefaultStream( ErrLog::Severity, const char*, int );

  /**
   *  Apply the logging test for a given severity level.
   *  The facility and code arguments are not used.
   */
  inline static bool _getDefaultLogging( ErrLog::Severity, const char*, int );

  /**
   *  Write a set of message parameters through into an ErrStream.
   *  This function exists so that subclasses of ErrLog can benefit
   *  from its friendship relationship with ErrStream.
   */
  static void _setParameters( ErrStream&, 
			      ErrLog::Severity, const char*, int );

  /**
   *  Evaluated and executes the default termination strategy for errors:
   *  abort() on fatal, exit(!=0) on severe, with the exit code derived from
   *  the message's "code" value.
   */
  static void _defaultTermination( ErrLog::Severity severity, 
				   int code );

  /**
   *  Tests the anticipated behavior of the default termination strategy.
   *  Available for use by implementations that need to know whether a
   *  termination will result before deciding on the disposition of a
   *  message.
   *
   *  @return 0 if termination will not occur, non-zero if it will.
   *          Negative if by signal (SIGABRT), positive if by exit(!=0).
   *
   *  @see QtrErrLog
   */
  static int _getDefaultTermination( ErrLog::Severity severity, 
				     int code );

//-----------------------------------------------------------------------
// Static private member functions --------------------------------------
//-----------------------------------------------------------------------

private:

  // These functions provide implementation details.

  static void _executeDefaultTermination( int value );


//-----------------------------------------------------------------------
// Static protected member data -----------------------------------------
//-----------------------------------------------------------------------

protected:

  // data members for the static front end and default implementation
  // These appear to be used in implementations.  A future cleanup pass
  // may be able to convert them all to private.
  static ErrStream* _myStreamCout;
  static ErrStream* _myStreamCerr;
  static ErrStream* _myDevnull;
  static Severity _minSeverity;


//-----------------------------------------------------------------------
// Static private member data -------------------------------------------
//-----------------------------------------------------------------------

private:

  // Singleton polymorphic implementation instance.
  static ErrLog* _implInstance;

  // data members for the static front end and default implementation
  static bool  _nestedWarnings;

private:

  bool _useStdio;

  // Contained by reference to avoid introducing header file dependencies:
  ErrOptRef<bool>* _useStdioOpt;
  ErrOptParser* _itsErrOptParser;

  // dis-allow copy construction and assignment
  ErrLog( const ErrLog& );
  ErrLog& operator=( const ErrLog& );

};


// Friends of ErrStream that are needed in implementing the "endmsg" feature.

void endmsg(ErrStream& es); // declare
std::ostream& operator<<(std::ostream& os, void (*fp)(ErrStream & es) );


// --------------------------------------------------------------------------
// Inline implementations ---------------------------------------------------
// --------------------------------------------------------------------------

inline bool 
ErrLog::logging( ErrLog::Severity severity, 
		 const char* facility, 
		 int code )
{
  return ( _implementation() != 0
	   ? _implementation()->doLogging( severity, facility, code )
	   : ErrLog::_getDefaultLogging( severity, facility, code )   );
}


inline void
ErrLog::_setImplementation( ErrLog* implementation )
{
  // "There can only be one" (or zero -- we can reset to the default).
  assert( implementation == 0 || _implementation() == 0 );

  // If the assertion is compiled away, this call is simply ignored.
  if ( implementation == 0 || _implementation() == 0 ) {
    _implInstance = implementation;
  }
}


inline bool 
ErrLog::_getDefaultLogging( ErrLog::Severity severity,
			    const char*,    // unused
			    int )           // unused
{
  return ( severity >= _minSeverity );
}

#endif // ERRLOG_HH
