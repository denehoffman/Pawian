//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universit??t Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//
#include <stdio.h>
#include <unistd.h>

#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <map>
#include <tuple>
#include <iterator>
#include <memory>


#include "PwaUtils/AbsLh.hh"
#include "AppUtils/AppBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/WelcomeScreen.hh"
#include "PwaUtils/DynRegistry.hh"

#include "ConfigParser/globalParser.hh"

#include "Utils/PawianCollectionUtils.hh"
#include "Utils/ErrLogUtils.hh"

#include "Event/EventReaderDefault.hh"
#include "Event/Event.hh"
#include "Event/EventList.hh"

#include "FitParams/ParamFactory.hh"
#include "FitParams/AbsPawianParameters.hh"

#include "ErrLogger/ErrLogger.hh"


int main(int __argc,char *__argv[]){
  char hostname[1024];
  gethostname(hostname, 1024);
  
  InfoMsg << welcomeScreen << endmsg;
  InfoMsg << "Compiled " << __DATE__ << " " << __TIME__ << endmsg;
  InfoMsg << "Hostname: " << hostname << endmsg;
  
  // Parse the command line
   globalParser* globalAppParams=new globalParser(__argc, __argv);

  // Set the desired error logging mode
  setErrLogMode(globalAppParams->getErrLogMode());

  // Setup the global environment and add the primary pbarp channel
  GlobalEnv::instance()->setup(globalAppParams);

  char* argvWoCfgFile[__argc];
  int argcWoCfgFile=0;
  for (int i=0; i<__argc ; ++i) {
    InfoMsg << "__argv[" << i << "]: " <<  __argv[i] << endmsg;
    std::string currentArgv(__argv[i]);
     if(currentArgv !=(char*)"--pbarpFiles" 
	&& currentArgv !=(char*)"--epemFiles"
	&& currentArgv !=(char*)"--ggFiles"
	&& currentArgv !=(char*)"--resFiles"
	&& currentArgv !=(char*)"--pipiScatteringFiles" 
	&& currentArgv !="-c" 
	&& currentArgv !="--configFile"){
      argvWoCfgFile[argcWoCfgFile]=__argv[i];
      argcWoCfgFile++;
    }
    else ++i;
  }

  AppBase theAppBase;
  theAppBase.addChannelEnvs(argcWoCfgFile, argvWoCfgFile);
  GlobalEnv::instance()->setupChannelEnvs(); 
 
  // Get mode
  std::string mode=globalAppParams->mode();
  
  theAppBase.createLhObjects();

  // Create likelihood objects
  ChannelEnvList channelEnvs=GlobalEnv::instance()->ChannelEnvs();

  //check replacements of parameter suffixes
  if (!GlobalEnv::instance()->areSuffixMapsIdentical()) return 0;
  
  //print out all replacements
  GlobalEnv::instance()->printFitParameterReplacements();
  
  if (mode=="dumpDefaultParams"){
    theAppBase.dumpDefaultParams();
    return 1;
  }
  
  if (mode=="dumpRandomParams"){
    theAppBase.dumpRandomParams();
    return 1;
  }

  // Read start param file
  std::shared_ptr<AbsPawianParameters> unsortedStartPawianParams=theAppBase.streamPawianParams();
  GlobalEnv::instance()->setStartPawianParams(unsortedStartPawianParams);
  std::shared_ptr<AbsPawianParameters> startPawianParams=GlobalEnv::instance()->startPawianParams();

  std::cout << "\n\n**************** Fit parameter **************************" << std::endl;
  startPawianParams->print(std::cout);

  // Fix params for all channels
  std::vector<std::string> fixedParams;

  //mnParFix must be set in global cfg file
//  std::vector<std::string> fixedGlobalParams = globalAppParams->fixedParams();
  std::vector<std::string> fixedGlobalParams = GlobalEnv::instance()->fixedParams();
  fixedParams.insert(fixedParams.end(), fixedGlobalParams.begin(), fixedGlobalParams.end());

  InfoMsg << "the fixed params are:" << endmsg; 
  std::vector<std::string>::iterator itStr;
  for (itStr=fixedParams.begin(); itStr!=fixedParams.end(); ++itStr) {
    InfoMsg << *itStr <<endmsg;
  }

  if(mode == "qaModeSimple")
    theAppBase.fixParams(startPawianParams,fixedParams, false);
  else
    theAppBase.fixParams(startPawianParams,fixedParams, true);

  //fill param list names for dynamics
  std::vector<std::shared_ptr<AbsDynamics> > dynVec=DynRegistry::instance()->getDynVec();
  std::vector<std::shared_ptr<AbsDynamics> >::iterator itDyn;
  for(itDyn=dynVec.begin(); itDyn!=dynVec.end(); ++itDyn) {
    (*itDyn)->fillParamNameList();
  }

  // Disable output buffering
  setvbuf(stdout, NULL, _IONBF, 0);

  if(mode == "client"){
    theAppBase.fitClientMode(startPawianParams);
  return 1;
 }

 if(mode == "server" || mode == "evoserver" || mode=="serverQA"){
   theAppBase.fitServerMode(startPawianParams);
   return 1;
 }

 if(mode == "qaModeSimple") {
   Alert << "qaModeSimple only working with singleChannelApp!!!" << endmsg;
   return 1;
 } else {
   Alert << "mode: " << mode << " is not supported in the coupled channel app!!!"<< endmsg;
   return 1;
 }

 return 1;
}


