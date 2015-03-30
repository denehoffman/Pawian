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

#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Particle/PdtParser.hh"

#include "PwaUtils/AbsLh.hh"
#include "FitParams/FitParColBase.hh"
#include "FitParams/StreamFitParColBase.hh"
#include "FitParams/PwaCovMatrix.hh"
#include "AppUtils/AppBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/NetworkClient.hh"
#include "PwaUtils/WelcomeScreen.hh"


#include "pbarpUtils/pbarpStatesLS.hh"
#include "ConfigParser/pbarpParser.hh"
#include "pbarpUtils/PbarpChannelEnv.hh"
#include "pbarpUtils/pbarpReaction.hh"
#include "pbarpUtils/pbarpBaseLh.hh"
#include "pbarpUtils/pbarpHeliLh.hh"
#include "pbarpUtils/pbarpCanoLh.hh"
#include "pbarpUtils/pbarpTensorLh.hh"
#include "pbarpUtils/pbarpHist.hh"
#include "pbarpUtils/spinDensityHist.hh"

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
  clock_t start, end;
  start= clock();

  for (int i=0; i<__argc ; ++i) Info << __argv[i] << endmsg;
  
  Info << welcomeScreen << endmsg;
  Info << "Compiled " << __DATE__ << " " << __TIME__ << endmsg;

  // Parse the command line
   globalParser* globalAppParams=new globalParser(__argc, __argv);

   std::vector<std::string> pbarpCfgs = globalAppParams->pbarpCfgs();
   std::vector<std::string> epemCfgs = globalAppParams->epemCfgs();

   //requirement single channel  sum reactionCfgs.size() == 1
   unsigned int numReactions=pbarpCfgs.size()+epemCfgs.size();
   Info << "numReactions: " << numReactions << endmsg;
   if (numReactions != 1){
     Alert << "for this single channel app it is required to define exactly 1 reaction!!!"
	   << "\n number of reactions here: " << numReactions << endmsg;
     exit(1);
   }

   GlobalEnv::instance()->setup(globalAppParams);
   ParserBase* channelParser=0; 
   std::shared_ptr<AbsChannelEnv> channelEnv;
   char* argvWoCfgFile[__argc];
   int argcWoCfgFile=0;
   for (int i=0; i<__argc ; ++i){
     Info << "__argv[" << i << "]: " <<  __argv[i] << endmsg;
     // if(__argv[i] ==(char*)"-c"){
     std::string currentArgv(__argv[i]);
     if(currentArgv =="-c" || currentArgv =="--configFile"){
       Alert << "for the singleCannelApp it is not allowed to use the flag -c !!!" << endmsg;
       exit(1); 
     }
     else if(currentArgv !=(char*)"--configPbarpFiles"){
       argvWoCfgFile[argcWoCfgFile]=__argv[i];
       argcWoCfgFile++;
     }
     else ++i;
   }

   bool isPbarpChannel=false;
   char* argvWCfgFile[argcWoCfgFile+2];
   for (int i=0; i<argcWoCfgFile ; ++i)  argvWCfgFile[i]=argvWoCfgFile[i];

   for(auto it=pbarpCfgs.begin(); it!=pbarpCfgs.end();++it){
    argvWCfgFile[argcWoCfgFile]=(char*)"-c";
    argvWCfgFile[argcWoCfgFile+1]=(char*)(*it).c_str();
    for (int i=0; i<argcWoCfgFile+2 ; ++i) Info << argvWCfgFile[i] << endmsg;
    pbarpParser* currentParser = new pbarpParser(argcWoCfgFile+2, argvWCfgFile);
     channelEnv = std::shared_ptr<AbsChannelEnv>(new PbarpChannelEnv(currentParser));
    channelParser=currentParser;
    GlobalEnv::instance()->AddEnv(channelEnv, AbsChannelEnv::CHANNEL_PBARP);
    isPbarpChannel=true;
  }
  for(auto it=epemCfgs.begin(); it!=epemCfgs.end();++it){
    argvWCfgFile[argcWoCfgFile]=(char*)"-c";
    argvWCfgFile[argcWoCfgFile+1]=(char*)(*it).c_str();
    epemParser* currentParser = new epemParser(argcWoCfgFile+2, argvWCfgFile);
    channelEnv = std::shared_ptr<AbsChannelEnv>(new EpemChannelEnv(currentParser));
    channelParser=currentParser;
    GlobalEnv::instance()->AddEnv(channelEnv, AbsChannelEnv::CHANNEL_EPEM);
  }

  GlobalEnv::instance()->replaceParser(channelParser);
  // Set the desired error logging mode
  setErrLogMode(channelParser->getErrLogMode());

  // Get mode
  std::string mode=channelParser->mode();
  
  AppBase theAppBase;
  theAppBase.createLhObjects();

  // Create likelihood objects
  //  ChannelEnvList channelEnvs=GlobalEnv::instance()->ChannelEnvs();

  if (mode=="dumpDefaultParams"){
    theAppBase.dumpDefaultParams();
    return 1;
  }

  // Read start param file
  std::string paramStreamerPath=channelParser->fitParamFile();
  std::string outputFileNameSuffix= GlobalEnv::instance()->outputFileNameSuffix();
  StreamFitParmsBase theParamStreamer(paramStreamerPath);
  fitParCol theStartparams=theParamStreamer.getFitParamVal();
  fitParCol theErrorparams=theParamStreamer.getFitParamErr();

  if (mode=="gen"){
    theAppBase.generate(theStartparams);
    return 1;
  }

  // Set minuit parameters
  std::shared_ptr<AbsPawianParameters> upar=ParamFactory::instance()->getParametersPointer("Minuit2");
  GlobalEnv::instance()->fitParColBase()->setAbsPawianParams(upar, theStartparams, theErrorparams);

  std::cout << "\n\n**************** Fit parameter **************************" << std::endl;
  for (int i=0; i<int(upar->Params().size()); ++i){
    std::cout << upar->Name(i) << "\t" << upar->Value(i) << "\t" << upar->Error(i) << std::endl;
  }

  // Fix params for all channels
  std::vector<std::string> fixedParams;
  std::vector<std::string> fixedChannelParams = channelParser->fixedParams();
  fixedParams.insert(fixedParams.end(), fixedChannelParams.begin(), fixedChannelParams.end());


  theAppBase.fixParams(upar,fixedParams);
  const unsigned int noOfFreeFitParams = upar->VariableParameters();

  // Disable output buffering
  setvbuf(stdout, NULL, _IONBF, 0);

  if(mode == "server" || mode == "evoserver"){
    theAppBase.fitServerMode(upar);
    return 1;
  }

  if(mode == "client"){

  std::ostringstream portStringStream;
  portStringStream << channelParser->serverPort();

  NetworkClient theClient(channelParser->serverAddress(), portStringStream.str());
  if(!theClient.Login())
    return 0;

  ChannelID channelID = theClient.channelID();

  bool cacheAmps = channelParser->cacheAmps();
  Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
  if (cacheAmps) GlobalEnv::instance()->Channel(channelID)->Lh()->cacheAmplitudes();

  const std::string datFile=channelParser->dataFile();
  const std::string mcFile=channelParser->mcFile();
  Info << "data file: " << datFile ;  // << endmsg;
  Info << "mc file: " << mcFile ;  // << endmsg;
  std::vector<std::string> dataFileNames;
  dataFileNames.push_back(datFile);

  std::vector<std::string> mcFileNames;
  mcFileNames.push_back(mcFile);

  EventList eventsDataClient;
  theAppBase.readEvents(eventsDataClient, dataFileNames, channelID, GlobalEnv::instance()->Channel(channelID)->useEvtWeight(), theClient.GetEventLimits()[0], theClient.GetEventLimits()[1]);

  EventList mcDataClient;
  theAppBase.readEvents(mcDataClient, mcFileNames, channelID, false, theClient.GetEventLimits()[2], theClient.GetEventLimits()[3]);

  std::shared_ptr<EvtDataBaseList> eventListPtr(new EvtDataBaseList(channelID));
  eventListPtr->read(eventsDataClient, mcDataClient);

  GlobalEnv::instance()->Channel(channelID)->Lh()->setDataVec(eventListPtr->getDataVecs());
  GlobalEnv::instance()->Channel(channelID)->Lh()->setMcVec(eventListPtr->getMcVecs());

  theAppBase.calcAndSendClientLh(theClient, theStartparams, channelID);

  return 1;
  }

  // The following modes only need the primary channel data/mc and lh ptr
  std::shared_ptr<AbsLh> theLhPtr = GlobalEnv::instance()->Channel()->Lh();

  const std::string datFile=channelParser->dataFile();
  const std::string mcFile=channelParser->mcFile();
  Info << "data file: " << datFile ;  // << endmsg;
  Info << "mc file: " << mcFile ;  // << endmsg;

  std::vector<std::string> dataFileNames;
  dataFileNames.push_back(datFile);

  std::vector<std::string> mcFileNames;
  mcFileNames.push_back(mcFile);

  if(mode == "spinDensity" && isPbarpChannel){
    bool cacheAmps = channelParser->cacheAmps();
    Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
    if (cacheAmps) GlobalEnv::instance()->Channel()->Lh()->cacheAmplitudes();
    
    EventList eventsData;
    theAppBase.readEvents(eventsData, dataFileNames, 0, GlobalEnv::instance()->Channel()->useEvtWeight(), 0, spinDensityHist::MAX_EVENTS);
    
    EventList mcData;
    theAppBase.readEvents(mcData, mcFileNames, 0, false, 0, spinDensityHist::MAX_EVENTS);
    
    std::shared_ptr<EvtDataBaseList> eventListPtr(new EvtDataBaseList(0));
    eventListPtr->read(eventsData, mcData);
    
    GlobalEnv::instance()->Channel()->Lh()->setDataVec(eventListPtr->getDataVecs());
    GlobalEnv::instance()->Channel()->Lh()->setMcVec(eventListPtr->getMcVecs());
    
    std::shared_ptr<spinDensityHist> theSpinDensityHist(new spinDensityHist(GlobalEnv::instance()->Channel()->Lh(), theStartparams));
    
    std::string serializationFileName = GlobalEnv::instance()->serializationFileName();
    std::ifstream serializationStream(serializationFileName.c_str());
    
    if(!serializationStream.is_open()){
      Warning << "Could not open serialization file." << endmsg;
    }
    else{
      boost::archive::text_iarchive boostInputArchive(serializationStream);
      std::shared_ptr<PwaCovMatrix> thePwaCovMatrix(new PwaCovMatrix);
      boostInputArchive >> *thePwaCovMatrix;
      theSpinDensityHist->SetCovarianceMatrix(thePwaCovMatrix);
    }
    theSpinDensityHist->Calculate();
    return 1;    
  }

  int noOfDataEvents = channelParser->noOfDataEvts();
  int ratioMcToData= channelParser->ratioMcToData();

  EventList eventsData;
  theAppBase.readEvents(eventsData, dataFileNames, 0, GlobalEnv::instance()->Channel()->useEvtWeight(), 0, noOfDataEvents);

  int maxMcEvts=eventsData.size()*ratioMcToData;

  EventList mcData;
  theAppBase.readEvents(mcData, mcFileNames, 0, false, 0, maxMcEvts-1);

  std::shared_ptr<EvtDataBaseList> eventListPtr(new EvtDataBaseList(0));

  if (mode=="plotMode"){
    theAppBase.plotMode(eventsData, mcData, eventListPtr);
    return 1;
  }

  if (mode=="qaModeSimple"){
    theAppBase.qaModeSimple(eventsData, mcData, theStartparams, eventListPtr, noOfFreeFitParams);
    return 1;
  }

  eventListPtr->read(eventsData, mcData);

  theLhPtr->setDataVec(eventListPtr->getDataVecs());
  theLhPtr->setMcVec(eventListPtr->getMcVecs());

  Info << "\nThe parameter values are: " << "\n" << endmsg;
  GlobalEnv::instance()->fitParColBase()->printParams(theStartparams);

  Info << "\nThe parameter errors are: " << "\n" << endmsg;
  GlobalEnv::instance()->fitParColBase()->printParams(theErrorparams);


  double evtWeightSumData = eventListPtr->NoOfWeightedDataEvts();
  if (mode=="qaMode"){
      theAppBase.qaMode(theStartparams, evtWeightSumData, noOfFreeFitParams );
      pbarpHist theHist;
      theHist.fillFromLhData(theLhPtr, theStartparams);
      end= clock();
      double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
      Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;

      return 1;
  }


  bool cacheAmps = channelParser->cacheAmps();
  Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
  if (cacheAmps) theLhPtr->cacheAmplitudes();

  if(mode=="pwa" || mode=="evo"){
    theAppBase.fitNonServerMode(upar, evtWeightSumData);
    return 1;
  }
  return 1;
}


