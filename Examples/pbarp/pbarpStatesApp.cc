#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include "ErrLogger/ErrLineLog.hh"
#include "PwaUtils/AbsStates.hh"
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

  std::vector< boost::shared_ptr<const jpcRes> > jpcStatesRequest;
  
  std::vector< boost::shared_ptr<const jpcRes> > theJPCStates=thepbarpState.jpcStates();
  ErrMsg(routine) << "The pbarp JPC states are: " << endmsg;
  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator it1;
  for ( it1=theJPCStates.begin(); it1!=theJPCStates.end(); ++it1){
    if (0!= (*it1))(*it1)->print(std::cout);
    std::cout << std::endl;

    if ( (*it1)->J==1 ||  (*it1)->J==4) jpcStatesRequest.push_back(*it1);
  }

  std::vector< boost::shared_ptr<const JPCLSM> >::const_iterator it;

  std::vector< boost::shared_ptr<const JPCLSM> > pbarpExtract=thepbarpState.extractJPCLSMStates(jpcStatesRequest);

  ErrMsg(routine) << "The requested pbarp states (JPCSML) are: " << endmsg;
  for ( it=pbarpExtract.begin(); it!=pbarpExtract.end(); ++it){
    if (0!= (*it))(*it)->print(std::cout);
  }

  std::vector< boost::shared_ptr<const JPCSM> >::const_iterator itJPCSM;

  std::vector< boost::shared_ptr<const JPCSM> > pbarpJPCSMExtract=thepbarpState.extractJPCSMStates(jpcStatesRequest);

  ErrMsg(routine) << "The requested pbarp states (JPCSM) are: " << endmsg;
  for ( itJPCSM=pbarpJPCSMExtract.begin(); itJPCSM!=pbarpJPCSMExtract.end(); ++itJPCSM){
    if (0!= (*itJPCSM))(*itJPCSM)->print(std::cout);
  }

  thepbarpState.print(std::cout);

  std::vector< boost::shared_ptr<const JPCLSM> > theSingletStates=thepbarpState.singletStates();
  ErrMsg(routine) << "The pbarp singlet states are: " << endmsg;
  for ( it=theSingletStates.begin(); it!=theSingletStates.end(); ++it){
    if (0!= (*it))(*it)->print(std::cout);
  }

  std::vector< boost::shared_ptr<const JPCLSM> > theTripletM0States=thepbarpState.tripletM0States();
  ErrMsg(routine) << "The pbarp triplet states with helicity=0 are: " << endmsg;
  for ( it=theTripletM0States.begin(); it!=theTripletM0States.end(); ++it){
    if (0!= (*it))(*it)->print(std::cout);
  }

  std::vector< boost::shared_ptr<const JPCLSM> > theTripletMp1States=thepbarpState.tripletMp1States();
  ErrMsg(routine) << "The pbarp triplet states with helicity=1 are: " << endmsg;
  for ( it=theTripletMp1States.begin(); it!=theTripletMp1States.end(); ++it){
    if (0!= (*it))(*it)->print(std::cout);
  }

  std::vector< boost::shared_ptr<const JPCLSM> > theTripletMm1States=thepbarpState.tripletMm1States();
  ErrMsg(routine) << "The pbarp triplet states with helicity=-1 are: " << endmsg;
  for ( it=theTripletMm1States.begin(); it!=theTripletMm1States.end(); ++it){
    if (0!= (*it))(*it)->print(std::cout);
  }

  if (0!=myLogger) delete myLogger;
  return 0;
}

