#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamParser.hh"


#include "ErrLogger/ErrLogger.hh"


int main(int __argc,char *__argv[]){

  // Parse the command line
  static Psi2STo2K2PiGamParser theAppParams(__argc, __argv);

  // Inform the audience about the execution mode
//   emitExecutionMode(theAppParams.getAppExecMode());

    
  return 0;
}

