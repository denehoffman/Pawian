#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include "ErrLogger/ErrLineLog.hh"
#include "PwaUtils/pbarpStates.hh"
#include "PwaUtils/DataUtils.hh"

int main(int __argc,char *__argv[]){

  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
		    || strcmp( __argv[1], "--help" ) == 0 ) ){

    std::cout << "test application which calculates all initial pbar p states\n"
	      << "output: JPC states with contributed LS combinations + helicicty + corresponding Clebsch-Gordan coefficient\n"
	      << "The application can be started with two flags: \n"
	      << "with -msg <errorLogMode>: choose the mode for the error logger \n"
	      << "with -jmax <value>: choose the maximum Spin J \n"
              << "e.g. type: ./pbarpStatesApp -jmax 6 -msg debugging \n"
	      << std::endl;
    return 0;
  }

  int optind=1;
  std::string msgModeStr="default";
  std::string jmaxStr="0";
  // decode arguments

  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-msg"){
      optind++;
      msgModeStr = __argv[optind];
      found=true;
    }
    if(sw=="-jmax"){
      optind++;
      jmaxStr = __argv[optind];
      found=true;
    }
    if (!found){
      ErrMsg(warning) << "Unknown switch: " 
            << __argv[optind] << endmsg;
      optind++;
    }
    
   while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
  }

  ErrLineLog* myLogger=0;  
  if(msgModeStr == "debugging") myLogger= new ErrLineLog(ErrLog::debugging);
  else if(msgModeStr == "trace") myLogger= new ErrLineLog(ErrLog::trace);
  else if(msgModeStr == "routine") myLogger= new ErrLineLog(ErrLog::routine);
  else if(msgModeStr == "warning")  myLogger= new ErrLineLog(ErrLog::warning);
  else if(msgModeStr == "error")    myLogger= new ErrLineLog(ErrLog::error); 
  else {
    myLogger= new ErrLineLog(ErrLog::routine);
    ErrMsg(warning) << "ErrorLogger not (properly) set -> Use mode 'ErrLog::routine' " << endmsg;  
  }



 std::stringstream jmaxStrStr(jmaxStr);
 int jmax=0;
 jmaxStrStr >> jmax;

  ErrMsg(routine) << "jmax: " << jmax << endmsg;

  pbarpStates thepbarpState(jmax);
  thepbarpState.print(std::cout);

  std::vector<PbarP*> theSingletStates=thepbarpState.singletStates();
  ErrMsg(routine) << "The pbarp singlet states are: " << endmsg;
  std::vector<PbarP*>::const_iterator it;
  for ( it=theSingletStates.begin(); it!=theSingletStates.end(); ++it){
    if (0!= (*it))(*it)->print(std::cout);
  }

  std::vector<PbarP*> theTripletM0States=thepbarpState.tripletM0States();
  ErrMsg(routine) << "The pbarp triplet states with helicity=0 are: " << endmsg;
  for ( it=theTripletM0States.begin(); it!=theTripletM0States.end(); ++it){
    if (0!= (*it))(*it)->print(std::cout);
  }

  std::vector<PbarP*> theTripletMp1States=thepbarpState.tripletMp1States();
  ErrMsg(routine) << "The pbarp triplet states with helicity=1 are: " << endmsg;
  for ( it=theTripletMp1States.begin(); it!=theTripletMp1States.end(); ++it){
    if (0!= (*it))(*it)->print(std::cout);
  }

  std::vector<PbarP*> theTripletMm1States=thepbarpState.tripletMm1States();
  ErrMsg(routine) << "The pbarp triplet states with helicity=-1 are: " << endmsg;
  for ( it=theTripletMm1States.begin(); it!=theTripletMm1States.end(); ++it){
    if (0!= (*it))(*it)->print(std::cout);
  }

  if (0!=myLogger) delete myLogger;
  return 0;
}

