#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKParser.hh"


#include "ErrLogger/ErrLogger.hh"


int main(int __argc,char *__argv[]){

  // Parse the command line
  static JpsiGamKsKlKKParser theAppParams(__argc, __argv);

  // Inform the audience about the execution mode
//   emitExecutionMode(theAppParams.getAppExecMode());

    
  return 0;
}

