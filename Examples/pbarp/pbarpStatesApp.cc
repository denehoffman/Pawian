#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/AbsStates.hh"
#include "pbarpUtils/pbarpStatesLS.hh"
#include "PwaUtils/DataUtils.hh"

int main(int __argc,char *__argv[]){

  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
		    || strcmp( __argv[1], "--help" ) == 0 ) ){

    std::cout << "test application which calculates all initial pbar p states\n"
	      << "output: JPC states with contributed LS combinations + helicicty + corresponding Clebsch-Gordan coefficient\n"
	      << "The application can be started with two flags: \n"
	      << "with -msg <errorLogMode>: choose the mode for the error logger \n"
	      << "with -lmax <value>: choose the maximum Spin J \n"
              << "e.g. type: ./pbarpStatesApp -lmax 6 -msg debugging \n"
	      << std::endl;
    return 0;
  }

  int optind=1;
  std::string msgModeStr="default";
  std::string lmaxStr="0";
  // decode arguments

  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-msg"){
      optind++;
      msgModeStr = __argv[optind];
      found=true;
    }
    if(sw=="-lmax"){
      optind++;
      lmaxStr = __argv[optind];
      found=true;
    }
    if (!found){
      Warning << "Unknown switch: " 
	      << __argv[optind] << endmsg;
      optind++;
    }
    
   while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
  }

  if(msgModeStr == "debugging") ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  else if(msgModeStr == "trace") ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
  else if(msgModeStr == "routine") ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
  else if(msgModeStr == "warning") ErrLogger::instance()->setLevel(log4cpp::Priority::WARN);
  else if(msgModeStr == "error")   ErrLogger::instance()->setLevel(log4cpp::Priority::ERROR);
  else {
    ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
    Warning << "ErrorLogger not (properly) set -> Use mode 'WARN' " << endmsg;  
  }



 std::stringstream lmaxStrStr(lmaxStr);
 int lmax=0;
 lmaxStrStr >> lmax;

  Info << "lmax: " << lmax << endmsg;

  pbarpStatesLS thepbarpState(lmax);

  thepbarpState.print(std::cout);

  return 0;
}

