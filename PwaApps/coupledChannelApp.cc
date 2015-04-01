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
#include "FitParams/FitParColBase.hh"
#include "FitParams/PwaCovMatrix.hh"
#include "AppUtils/AppBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/NetworkClient.hh"
#include "PwaUtils/WelcomeScreen.hh"


#include "pbarpUtils/pbarpStatesLS.hh"
#include "ConfigParser/pbarpParser.hh"
#include "pbarpUtils/PbarpChannelEnv.hh"

#include "ConfigParser/globalParser.hh"
#include "ConfigParser/epemParser.hh"
#include "epemUtils/EpemChannelEnv.hh"

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

  // Create environments for coupled channels
  // Currently secondary channels need to be of the same type(pbarp, epem,..) as the primary one
  std::vector<std::string> pbarpCfgs = globalAppParams->pbarpCfgs();
  for(auto it=pbarpCfgs.begin(); it!=pbarpCfgs.end();++it){
    char* argv[] = {__argv[0], (char*)"-c", (char*)(*it).c_str()};
    pbarpParser* ccParser = new pbarpParser(3, argv);
    GlobalEnv::instance()->AddEnv(std::shared_ptr<PbarpChannelEnv>(new PbarpChannelEnv(ccParser)), AbsChannelEnv::CHANNEL_PBARP);
  }

  std::vector<std::string> epemCfgs = globalAppParams->epemCfgs();
  for(auto it=epemCfgs.begin(); it!=epemCfgs.end();++it){
    char* argv[] = {__argv[0], (char*)"-c", (char*)(*it).c_str()};
    epemParser* ccParser = new epemParser(3, argv);
    GlobalEnv::instance()->AddEnv(std::shared_ptr<EpemChannelEnv>(new EpemChannelEnv(ccParser)), AbsChannelEnv::CHANNEL_EPEM);
  }

  // Get mode
  std::string mode=globalAppParams->mode();
  
  AppBase theAppBase;
  theAppBase.createLhObjects();

  // Create likelihood objects
  ChannelEnvList channelEnvs=GlobalEnv::instance()->ChannelEnvs();

  if (mode=="dumpDefaultParams"){
    theAppBase.dumpDefaultParams();
    return 1;
  }

  // Read start param file
  fitParCol theStartparams;
  fitParCol theErrorparams;
  theAppBase.streamParams(theStartparams, theErrorparams); 

  // Set minuit parameters
  std::shared_ptr<AbsPawianParameters> upar=ParamFactory::instance()->getParametersPointer("Minuit2");
  GlobalEnv::instance()->fitParColBase()->setAbsPawianParams(upar, theStartparams, theErrorparams);

  std::cout << "\n\n**************** Fit parameter **************************" << std::endl;
  for (int i=0; i<int(upar->Params().size()); ++i){
    std::cout << upar->Name(i) << "\t" << upar->Value(i) << "\t" << upar->Error(i) << std::endl;
  }

  // Fix params for all channels
  std::vector<std::string> fixedParams;
  for(auto it=channelEnvs.begin();it!=channelEnvs.end();++it){
    std::vector<std::string> fixedChannelParams = (*it).first->parser()->fixedParams();
    fixedParams.insert(fixedParams.end(), fixedChannelParams.begin(), fixedChannelParams.end());
  }
  theAppBase.fixParams(upar,fixedParams);

  // Disable output buffering
  setvbuf(stdout, NULL, _IONBF, 0);

  if(mode == "client"){
    theAppBase.fitClientMode(theStartparams);
  return 1;
 }

 if(mode == "server" || mode == "evoserver"){
   theAppBase.fitServerMode(upar);
   return 1;
 }

 return 1;
}


