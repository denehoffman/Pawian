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

  Info << welcomeScreen << endmsg;
  Info << "Compiled " << __DATE__ << " " << __TIME__ << endmsg;

  // Parse the command line
   globalParser* globalAppParams=new globalParser(__argc, __argv);

  // Set the desired error logging mode
  setErrLogMode(globalAppParams->getErrLogMode());

  // Setup the global environment and add the primary pbarp channel
  GlobalEnv::instance()->setup(globalAppParams);

  char* argv[3];
  int argc=3;
  AppBase theAppBase;
  theAppBase.addChannelEnvs(argc, argv);
  GlobalEnv::instance()->setupChannelEnvs(); 
 
  // Get mode
  std::string mode=globalAppParams->mode();
  
  theAppBase.createLhObjects();

  // Create likelihood objects
  ChannelEnvList channelEnvs=GlobalEnv::instance()->ChannelEnvs();

  if (mode=="dumpDefaultParams"){
    theAppBase.dumpDefaultParams();
    return 1;
  }

  // Read start param file
  std::shared_ptr<AbsPawianParameters> startPawianParams=theAppBase.streamPawianParams();
  // fitParCol theStartparams;
  // fitParCol theErrorparams;
  // theAppBase.streamParams(theStartparams, theErrorparams); 

  // Set minuit parameters
  // std::shared_ptr<AbsPawianParameters> upar=ParamFactory::instance()->getParametersPointer("Minuit2");
  // GlobalEnv::instance()->fitParColBase()->setAbsPawianParams(upar, theStartparams, theErrorparams);

  std::cout << "\n\n**************** Fit parameter **************************" << std::endl;
  startPawianParams->print(std::cout);
  // for (int i=0; i<int(upar->Params().size()); ++i){
  //   std::cout << upar->Name(i) << "\t" << upar->Value(i) << "\t" << upar->Error(i) << std::endl;
  // }

  // Fix params for all channels
  std::vector<std::string> fixedParams;
  for(auto it=channelEnvs.begin();it!=channelEnvs.end();++it){
    std::vector<std::string> fixedChannelParams = (*it).first->parser()->fixedParams();
    fixedParams.insert(fixedParams.end(), fixedChannelParams.begin(), fixedChannelParams.end());
  }
  theAppBase.fixParams(startPawianParams,fixedParams);

  // Disable output buffering
  setvbuf(stdout, NULL, _IONBF, 0);

  if(mode == "client"){
    theAppBase.fitClientMode(startPawianParams);
  return 1;
 }

 if(mode == "server" || mode == "evoserver"){
   theAppBase.fitServerMode(startPawianParams);
   return 1;
 }

 return 1;
}


