#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiParser.hh"


#include "ErrLogger/ErrLogger.hh"


int main(int __argc,char *__argv[]){

  // Parse the command line
  static JpsiGamEtaPiPiParser theAppParams(__argc, __argv);

  // Inform the audience about the execution mode
//   emitExecutionMode(theAppParams.getAppExecMode());

    
  return 0;
}

