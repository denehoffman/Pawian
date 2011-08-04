#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamParser.hh"


#include "ErrLogger/ErrLogger.hh"


int main(int __argc,char *__argv[]){

  // Parse the command line
  static Psi2SToKpKmPiGamParser theAppParams(__argc, __argv);

  // Inform the audience about the execution mode
//   emitExecutionMode(theAppParams.getAppExecMode());

    
  return 0;
}

