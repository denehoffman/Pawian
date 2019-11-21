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
#include "FitParams/PwaCovMatrix.hh"
#include "AppUtils/AppBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/EvtDataListFactory.hh"
#include "PwaUtils/WelcomeScreen.hh"
#include "PwaUtils/DynRegistry.hh"

#include "pbarpUtils/spinDensityHist.hh"

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
  clock_t start, end;
  start= clock();

  for (int i=0; i<__argc ; ++i) {
    InfoMsg << __argv[i] << endmsg;
  }

  char hostname[1024];
  gethostname(hostname, 1024); 
  InfoMsg << welcomeScreen << endmsg;
  InfoMsg << "Compiled " << __DATE__ << " " << __TIME__ << endmsg;
  InfoMsg << "Hostname: " << hostname << endmsg;

  // Parse the command line
  globalParser* globalAppParams=new globalParser(__argc, __argv);

  std::vector<std::string> pbarpCfgs = globalAppParams->pbarpCfgs();
  std::vector<std::string> epemCfgs = globalAppParams->epemCfgs();
  std::vector<std::string> resCfgs = globalAppParams->resCfgs();
  std::vector<std::string> ggCfgs = globalAppParams->ggCfgs();
  std::vector<std::string> pipiScatteringCfgs = globalAppParams->pipiScatteringCfgs();
  
  //requirement single channel  sum reactionCfgs.size() == 1
  unsigned int numReactions=pbarpCfgs.size()+epemCfgs.size()+resCfgs.size()+ggCfgs.size()+pipiScatteringCfgs.size();
  InfoMsg << "numReactions: " << numReactions << endmsg;
  if (numReactions != 1){
    Alert << "for this single channel app it is required to define exactly 1 reaction!!!"
	  << "\n number of reactions here: " << numReactions << endmsg;
    exit(1);
  }
  
  GlobalEnv::instance()->setup(globalAppParams);
  
  char* argvWoCfgFile[__argc];
  int argcWoCfgFile=0;
  for (int i=0; i<__argc ; ++i) {
    InfoMsg << "__argv[" << i << "]: " <<  __argv[i] << endmsg;
    std::string currentArgv(__argv[i]);
    if(currentArgv =="-c" || currentArgv =="--configFile"){
      Alert << "for the singleCannelApp it is not allowed to use the flag -c !!!" << endmsg;
      exit(1); 
    }
    else if(currentArgv !=(char*)"--pbarpFiles" && currentArgv !=(char*)"--epemFiles" 
	    && currentArgv !=(char*)"--resFiles" && currentArgv !=(char*)"--ggFiles" 
	    && currentArgv !=(char*)"--pipiScatteringFiles"){
      argvWoCfgFile[argcWoCfgFile]=__argv[i];
      argcWoCfgFile++;
    }
    else ++i;
  }
  
  bool isPbarpChannel=false;
  bool isPiPiScatteringChannel=false;
  if (pbarpCfgs.size()==1) isPbarpChannel=true;
  else if (pipiScatteringCfgs.size()==1) isPiPiScatteringChannel=true;
  
  AppBase theAppBase;
  theAppBase.addChannelEnvs(argcWoCfgFile, argvWoCfgFile);
  
  GlobalEnv::instance()->replaceParser(GlobalEnv::instance()->Channel(0)->parser());
  std::vector<int> channelIDs; 
  channelIDs.push_back(GlobalEnv::instance()->parser()->singleChannelId());
  
  Alert << "trying to set channel ID to " << channelIDs[0] << endmsg;
  GlobalEnv::instance()->setupChannelEnvs(channelIDs);
  
  // Set the desired error logging mode
  setErrLogMode(GlobalEnv::instance()->parser()->getErrLogMode());
  
  // Get mode
  std::string mode=GlobalEnv::instance()->parser()->mode();
  
  theAppBase.createLhObjects();
  
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
  
  
  if (mode=="gen"){
    if(isPiPiScatteringChannel){
      Alert << "gen mode is not supported for pipi scattering reactions!!!" << endmsg;
      exit(1);
    }

    theAppBase.generate(startPawianParams);
    return 1;
  }

  std::cout << "\n\n**************** Fit parameter **************************" << std::endl;
  startPawianParams->print(std::cout);
  // for (int i=0; i<int(upar->Params().size()); ++i){
  //   std::cout << upar->Name(i) << "\t" << upar->Value(i) << "\t" << upar->Error(i) << std::endl;
  // }

  // Fix params for all channels
  if (GlobalEnv::instance()->parser()->doScaling()){
    if(isPiPiScatteringChannel){
      Alert << "scaling fits are not supported for pipi scattering reactions!!!" << endmsg;
      exit(1);
    }

    theAppBase.fixAllReleaseScaleParams(startPawianParams);
  }
  else{
    std::vector<std::string> fixedParams;
    std::vector<std::string> fixedChannelParams = GlobalEnv::instance()->parser()->fixedParams();
    fixedParams.insert(fixedParams.end(), fixedChannelParams.begin(), fixedChannelParams.end());
    theAppBase.fixParams(startPawianParams,fixedParams);
  }

  //fill param list names for dynamics
  std::vector<std::shared_ptr<AbsDynamics> > dynVec=DynRegistry::instance()->getDynVec();
  std::vector<std::shared_ptr<AbsDynamics> >::iterator itDyn;
  for(itDyn=dynVec.begin(); itDyn!=dynVec.end(); ++itDyn){
    (*itDyn)->fillParamNameList();
  } 

  // Disable output buffering
  setvbuf(stdout, NULL, _IONBF, 0);

  if(mode == "server" || mode == "evoserver"){
    theAppBase.fitServerMode(startPawianParams);
    return 1;
  }

  if(mode == "client"){
    theAppBase.fitClientMode(startPawianParams);
  return 1;
  }


  // The following modes only need the primary channel data/mc and lh ptr
  std::shared_ptr<AbsLh> theLhPtr = GlobalEnv::instance()->Channel()->Lh();

  const std::string datFile=GlobalEnv::instance()->parser()->dataFile();
  InfoMsg << "data file: " << datFile << endmsg;
  std::vector<std::string> dataFileNames;
  dataFileNames.push_back(datFile);


  const std::string mcFile=GlobalEnv::instance()->parser()->mcFile();
  std::vector<std::string> mcFileNames;
  if(!isPiPiScatteringChannel) {
    InfoMsg << "mc file: " << mcFile << endmsg;
    mcFileNames.push_back(mcFile);
  }

  std::shared_ptr<EvtDataBaseList> eventListPtr = 
    EvtDataListFactory::instance()->evtDataListPtr(GlobalEnv::instance()->Channel());
  
  if(mode == "spinDensity" && isPbarpChannel){
    bool cacheAmps = GlobalEnv::instance()->parser()->cacheAmps();
    InfoMsg << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
    if (cacheAmps) GlobalEnv::instance()->Channel()->Lh()->cacheAmplitudes();
    
    EventList eventsData;
    theAppBase.readEvents(eventsData, dataFileNames, 0, GlobalEnv::instance()->Channel()->useDataEvtWeight(), 
			  0, spinDensityHist::MAX_EVENTS);
    
    EventList mcData;
    theAppBase.readEvents(mcData, mcFileNames, 0, GlobalEnv::instance()->Channel()->useMCEvtWeight(), 0, spinDensityHist::MAX_EVENTS);
    
    eventListPtr->read(eventsData, mcData);
    
    GlobalEnv::instance()->Channel()->Lh()->setDataVec(eventListPtr->getDataVecs());
    GlobalEnv::instance()->Channel()->Lh()->setMcVec(eventListPtr->getMcVecs());
    
    std::shared_ptr<spinDensityHist> theSpinDensityHist(new spinDensityHist(GlobalEnv::instance()->Channel()->Lh(), startPawianParams));
    
    std::string serializationFileName = GlobalEnv::instance()->serializationFileName();
    std::ifstream serializationStream(serializationFileName.c_str());
    
    if(!serializationStream.is_open()) {
      WarningMsg << "Could not open serialization file." << endmsg;
    } else {
      boost::archive::text_iarchive boostInputArchive(serializationStream);
      std::shared_ptr<PwaCovMatrix> thePwaCovMatrix(new PwaCovMatrix);
      boostInputArchive >> *thePwaCovMatrix;
      theSpinDensityHist->SetCovarianceMatrix(thePwaCovMatrix);
    }
    theSpinDensityHist->Calculate();
    return 1;    
  }

  int noOfDataEvents = GlobalEnv::instance()->parser()->noOfDataEvts();
  int ratioMcToData= GlobalEnv::instance()->parser()->ratioMcToData();

  EventList eventsData;
  EventList mcData;
  if(isPiPiScatteringChannel) {
    theAppBase.readScatteringEvents(eventsData, dataFileNames, 0);
  } else {
    theAppBase.readEvents(eventsData, dataFileNames, 0, GlobalEnv::instance()->Channel()->useDataEvtWeight(), 
			  0, noOfDataEvents);
    
    int maxMcEvts=eventsData.size()*ratioMcToData;

    theAppBase.readEvents(mcData, mcFileNames, 0, GlobalEnv::instance()->Channel()->useMCEvtWeight(), 0, maxMcEvts-1);
  }

  if (mode=="plotMode"){
    theAppBase.plotMode(eventsData, mcData, eventListPtr);
    return 1;
  }

  if (mode=="qaModeSimple"){
    if(isPiPiScatteringChannel){
      Alert << "mode qaModeSimple is not supported for pipi scattering reactions!!!" 
	    << "\nuse qaMode instead!!!" << endmsg;
      exit(1);
    }
    theAppBase.qaModeSimple(eventsData, mcData, startPawianParams);
    return 1;
  }

  if (mode=="qaModeEffCorrection"){
    if(isPiPiScatteringChannel){
      Alert << "mode qaModeEffCorrection is not supported for pipi scattering reactions!!!" 
	    << endmsg;
      exit(1);
    }
    const std::string truthFile=GlobalEnv::instance()->parser()->truthFile();
    int ratioTruthToMc= GlobalEnv::instance()->parser()->ratioTruthToMc();
    int maxTruthEvts=eventsData.size()*ratioMcToData*ratioTruthToMc;
    InfoMsg << "truth file: " << truthFile << endmsg;
    std::vector<std::string> truthFileNames;
    truthFileNames.push_back(truthFile);
    EventList truthData;
    theAppBase.readEvents(truthData, truthFileNames, 0, GlobalEnv::instance()->parser()->useTruthEvtWeight(), 
			  0, maxTruthEvts-1);
    theAppBase.qaModeEffCorrection(eventsData, mcData, truthData, startPawianParams);
    return 1;
  }

  eventListPtr->read(eventsData, mcData);

  eventsData.removeAndDeleteEvents(0, eventsData.size()-1);
  mcData.removeAndDeleteEvents(0, mcData.size()-1);

  theLhPtr->setDataVec(eventListPtr->getDataVecs());
  theLhPtr->setMcVec(eventListPtr->getMcVecs());

  InfoMsg << "\nThe parameter values and errors are: " << "\n" << endmsg;
  startPawianParams->print(std::cout);

  double evtWeightSumData = eventListPtr->NoOfWeightedDataEvts();
  InfoMsg << "evtWeightSumData: " << evtWeightSumData << endmsg;
  double evtWeightSumMc = eventListPtr->NoOfWeightedMcEvts();
  InfoMsg << "evtWeightSumMc: " << evtWeightSumMc << endmsg;

  if (mode=="qaMode"){
      theAppBase.qaMode(startPawianParams, evtWeightSumData );
      end= clock();
      double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
      InfoMsg << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;

      return 1;
  }

  bool cacheAmps = GlobalEnv::instance()->parser()->cacheAmps();
  InfoMsg << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
  if (cacheAmps) theLhPtr->cacheAmplitudes();

  if(mode=="pwa" || mode=="evo"){
    theAppBase.fitNonServerMode(startPawianParams, evtWeightSumData, evtWeightSumMc);
    return 1;
  }
  return 1;
}


