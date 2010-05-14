//--------------------------------------------------------------------------
//
// Description:
//
//  Error & logging with line-number class test routine
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Bob Jacobsen
//      $Id: ErrLineLogTest.cc,v 1.1.1.1 2005/03/29 17:04:19 steinke Exp $
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------

//-----------------------
// Test Class's Header --
//-----------------------
#include "ErrLogger/ErrLineLog.hh"

#include <iostream>
using std::cout;
using std::endl;
using std::ostream;

class testOp {
private:
  friend ostream& operator<<(ostream & os,const testOp& tO) {
    os << "friend operator<< called";  return os;
  };
}; 

int main() {
    // construct an ErrLogger with a minimum severity
    ErrLineLog myLogger(ErrLog::warning);

    cout << "First, log with literal facility string, code integer" << endl;
    ErrLog::msg(ErrLog::warning, "someFacility", 703 )<<"some text"<<endmsg;

    cout << "Now with code 0" << endl;
    ErrLog::msg(ErrLog::warning, "someFacility", 0 )<<"some text"<<endmsg;

    cout << "Now log through string-code interface" << endl;
    ErrLog::msg(ErrLog::warning, "someFacility", "code" )<<"some text"<<endmsg;

    cout << "Now try using ErrMsg macro"<< endl;
    ErrMsg(warning) << "tracing through macro"<< endmsg;

    cout << "Now see what happens with nested messages" << endl;
    ErrMsg(error) << "incomplete message at error severity" << endl;
    ErrMsg(routine) << "subsequent message at routine severity" << endmsg;

    testOp myOp;
    ErrMsg(warning) << myOp << endmsg;
 
    return 0;
}
