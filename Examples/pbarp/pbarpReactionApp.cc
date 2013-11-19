//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
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

#include "TROOT.h"

#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Particle/PdtParser.hh"

#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/StreamFitParmsBase.hh"
#include "PwaUtils/PwaFcnBase.hh"
#include "PwaUtils/PwaFcnServer.hh"
#include "PwaUtils/PwaCovMatrix.hh"
#include "PwaUtils/WaveContribution.hh"
//#include "PwaUtils/PwaGen.hh"
#include "PwaUtils/AppBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/EvtWeightList.hh"
#include "PwaUtils/NetworkServer.hh"
#include "PwaUtils/NetworkClient.hh"
#include "PwaUtils/WelcomeScreen.hh"
#include "PwaUtils/EvoMinimizer.hh"

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

#include "Utils/PawianCollectionUtils.hh"
#include "Utils/ErrLogUtils.hh"

#include "Event/EventReaderDefault.hh"
#include "Event/Event.hh"
#include "Event/EventList.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnScan.h"

#include "ErrLogger/ErrLogger.hh"



int main(int __argc,char *__argv[]){
  clock_t start, end;
  start= clock();

  Info << welcomeScreen << endmsg;
  Info << "Compiled " << __DATE__ << " " << __TIME__ << endmsg;

  // Parse the command line
  pbarpParser* theAppParams=new pbarpParser(__argc, __argv);

  // Set the desired error logging mode
  setErrLogMode(theAppParams->getErrLogMode());

  // Setup the global environment and add the primary pbarp channel
  GlobalEnv::instance()->setup(theAppParams);
  GlobalEnv::instance()->AddEnv(std::shared_ptr<PbarpChannelEnv>(new PbarpChannelEnv(theAppParams)), AbsChannelEnv::CHANNEL_PBARP);

  // Print the primary pbarpReaction
  std::shared_ptr<pbarpReaction> thepbarpReaction=std::static_pointer_cast<PbarpChannelEnv>(GlobalEnv::instance()->PbarpChannel())->reaction();
  thepbarpReaction->print(std::cout);

  // Create environments for coupled channels
  // Currently secondary channels need to be of the same type(pbarp, epem,..) as the primary one
  std::vector<std::string> coupledChannelCfgs = theAppParams->coupledChannelCfgs();
  for(auto it=coupledChannelCfgs.begin(); it!=coupledChannelCfgs.end();++it){
     char* argv[] = {__argv[0], (char*)"-c", (char*)(*it).c_str()};
     pbarpParser* ccParser = new pbarpParser(3, argv);
     GlobalEnv::instance()->AddEnv(std::shared_ptr<PbarpChannelEnv>(new PbarpChannelEnv(ccParser)), AbsChannelEnv::CHANNEL_PBARP);
  }

  // Get mode
  std::string mode=theAppParams->mode();

  // Create likelihood objects
  ChannelEnvList channelEnvs=GlobalEnv::instance()->ChannelEnvs();
  for(auto it=channelEnvs.begin();it!=channelEnvs.end();++it){
     std::string prodFormalism=(*it).first->parser()->productionFormalism();
     std::shared_ptr<AbsLh> theLhPtr;
     if(prodFormalism=="Cano") theLhPtr=std::shared_ptr<AbsLh>(new pbarpCanoLh((*it).first->channelID()));
     else if(prodFormalism=="Heli") theLhPtr=std::shared_ptr<AbsLh>(new pbarpHeliLh((*it).first->channelID()));
     else if(prodFormalism=="Tensor") theLhPtr=std::shared_ptr<AbsLh>(new pbarpTensorLh((*it).first->channelID()));
     else {
         Alert << "prodFormalism\t" << prodFormalism << "\tdoesn't exist!!!" << endmsg;
         exit(1);
     }
     (*it).first->SetLh(theLhPtr);
  }

  // Generate the full parameter set using the likelihood list
  GlobalEnv::instance()->CreateDefaultParameterSet();




  AppBase theAppBase;

  if (mode=="dumpDefaultParams"){
    theAppBase.dumpDefaultParams();
    return 1;
  }

  // Read start param file
  std::string paramStreamerPath=theAppParams->fitParamFile();
  std::string outputFileNameSuffix= GlobalEnv::instance()->outputFileNameSuffix();
  StreamFitParmsBase theParamStreamer(paramStreamerPath);
  fitParams theStartparams=theParamStreamer.getFitParamVal();
  fitParams theErrorparams=theParamStreamer.getFitParamErr();

  if (mode=="gen"){
    theAppBase.generate(theStartparams);
    return 1;
  }

  // Set minuit parameters
  MnUserParameters upar;
  GlobalEnv::instance()->fitParamsBase()->setMnUsrParams(upar, theStartparams, theErrorparams);

  std::cout << "\n\n**************** Minuit Fit parameter **************************" << std::endl;
  for (int i=0; i<int(upar.Params().size()); ++i){
    std::cout << upar.Name(i) << "\t" << upar.Value(i) << "\t" << upar.Error(i) << std::endl;
  }

  // Fix params for all channels
  std::vector<std::string> fixedParams;
  for(auto it=channelEnvs.begin();it!=channelEnvs.end();++it){
      std::vector<std::string> fixedChannelParams = (*it).first->parser()->fixedParams();
      fixedParams.insert(fixedParams.end(), fixedChannelParams.begin(), fixedChannelParams.end());
  }
  theAppBase.fixParams(upar,fixedParams);
  const unsigned int noOfFreeFitParams = upar.VariableParameters();

  // Disable output buffering
  setvbuf(stdout, NULL, _IONBF, 0);

  if(mode == "client"){

  std::ostringstream portStringStream;
  portStringStream << theAppParams->serverPort();

  NetworkClient theClient(theAppParams->serverAddress(), portStringStream.str());
  if(!theClient.Login())
    return 0;

  ChannelID channelID = theClient.channelID();

  bool cacheAmps = theAppParams->cacheAmps();
  Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
  if (cacheAmps) GlobalEnv::instance()->Channel(channelID)->Lh()->cacheAmplitudes();

  const std::string datFile=GlobalEnv::instance()->Channel(channelID)->parser()->dataFile();
  const std::string mcFile=GlobalEnv::instance()->Channel(channelID)->parser()->mcFile();
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

  std::shared_ptr<EvtDataBaseList> pbarpEventListPtr(new EvtDataBaseList(channelID));
  pbarpEventListPtr->read(eventsDataClient, mcDataClient);

  GlobalEnv::instance()->Channel(channelID)->Lh()->setDataVec(pbarpEventListPtr->getDataVecs());
  GlobalEnv::instance()->Channel(channelID)->Lh()->setMcVec(pbarpEventListPtr->getMcVecs());

  theAppBase.calcAndSendClientLh(theClient, theStartparams, channelID);

  return 1;
 }



  // Get the mc data limit
  int ratioMcToData=theAppParams->ratioMcToData();



 if(mode == "server"){

    double evtWeightSumData=0;
    ChannelEnvList channelEnvs=GlobalEnv::instance()->ChannelEnvs();
    std::map<short, std::tuple<long, double, long> > numEventMap;

    for(auto it=channelEnvs.begin();it!=channelEnvs.end();++it){
      const std::string datFile=(*it).first->parser()->dataFile();
      const std::string mcFile=(*it).first->parser()->mcFile();
      Info << "data file: " << datFile ;  // << endmsg;
      Info << "mc file: " << mcFile ;  // << endmsg;
      std::vector<std::string> dataFileNames;
      dataFileNames.push_back(datFile);

      std::vector<std::string> mcFileNames;
      mcFileNames.push_back(mcFile);

      EventList eventsData;
      theAppBase.readEvents(eventsData, dataFileNames, (*it).first->channelID(), (*it).first->useEvtWeight());

      EventList mcData;
      int maxMcEvts=eventsData.size()*ratioMcToData;
      theAppBase.readEvents(mcData, mcFileNames, (*it).first->channelID(), false, 0, maxMcEvts-1);

      std::shared_ptr<EvtWeightList> pbarpWeightListPtr(new EvtWeightList((*it).first->channelID()));
      pbarpWeightListPtr->read(eventsData, mcData);
      evtWeightSumData+=pbarpWeightListPtr->NoOfWeightedDataEvts();

      numEventMap[(*it).first->channelID()] = std::tuple<long, double,long>(eventsData.size(), pbarpWeightListPtr->NoOfWeightedDataEvts(), mcData.size());
    }

    std::shared_ptr<NetworkServer> theServer(new NetworkServer(theAppParams->serverPort(), theAppParams->noOfClients(), numEventMap));

    PwaFcnServer theFcnServer(theServer);
    theServer->WaitForFirstClientLogin();

    FunctionMinimum min=theAppBase.migradDefault(theFcnServer, upar);

    theServer->BroadcastClosingMessage();
    Info << "Closing server." << endmsg;

    theAppBase.printFitResult(min, theStartparams, std::cout, evtWeightSumData, noOfFreeFitParams);

    return 1;
 }


 if(mode == "evoserver"){
    double evtWeightSumData=0;
    ChannelEnvList channelEnvs=GlobalEnv::instance()->ChannelEnvs();
    std::map<short, std::tuple<long, double, long> > numEventMap;

    for(auto it=channelEnvs.begin();it!=channelEnvs.end();++it){
      const std::string datFile=(*it).first->parser()->dataFile();
      const std::string mcFile=(*it).first->parser()->mcFile();
      Info << "data file: " << datFile ;  // << endmsg;
      Info << "mc file: " << mcFile ;  // << endmsg;
      std::vector<std::string> dataFileNames;
      dataFileNames.push_back(datFile);

      std::vector<std::string> mcFileNames;
      mcFileNames.push_back(mcFile);

      EventList eventsData;
      theAppBase.readEvents(eventsData, dataFileNames, (*it).first->channelID(), (*it).first->useEvtWeight());

      EventList mcData;
      int maxMcEvts=eventsData.size()*ratioMcToData;
      theAppBase.readEvents(mcData, mcFileNames, (*it).first->channelID(), false, 0, maxMcEvts-1);

      std::shared_ptr<EvtWeightList> pbarpWeightListPtr(new EvtWeightList((*it).first->channelID()));
      pbarpWeightListPtr->read(eventsData, mcData);
      evtWeightSumData+=pbarpWeightListPtr->NoOfWeightedDataEvts();

      numEventMap[(*it).first->channelID()] = std::tuple<long, double,long>(eventsData.size(), pbarpWeightListPtr->NoOfWeightedDataEvts(), mcData.size());
    }

   std::shared_ptr<NetworkServer> theServer(new NetworkServer(theAppParams->serverPort(), theAppParams->noOfClients(), numEventMap));


   PwaFcnServer theFcnServer(theServer);
   theServer->WaitForFirstClientLogin();

   EvoMinimizer theEvoMinimizer(theFcnServer, upar, GlobalEnv::instance()->parser()->evoPopulation(), GlobalEnv::instance()->parser()->evoIterations());
   Info <<"start evolutionary minimizer "<< endmsg;
   std::vector<double> finalParamVec = theEvoMinimizer.Minimize();

   theServer->BroadcastClosingMessage();
   Info << "Closing server." << endmsg;

   fitParams finalFitParams=theStartparams;
   GlobalEnv::instance()->fitParamsBase()->getFitParamVal(finalParamVec, finalFitParams);

   fitParams finalFitErrs=theErrorparams;

   std::ostringstream finalResultname;
   finalResultname << "finalResult" << outputFileNameSuffix << ".dat";

   std::ofstream theStream ( finalResultname.str().c_str() );
   GlobalEnv::instance()->fitParamsBase()->dumpParams(theStream, finalFitParams, finalFitErrs);

   return 1;
 }






  // The following modes only need the primary channel data/mc and lh ptr
  std::shared_ptr<AbsLh> theLhPtr = GlobalEnv::instance()->Channel()->Lh();

  const std::string datFile=theAppParams->dataFile();
  const std::string mcFile=theAppParams->mcFile();
  Info << "data file: " << datFile ;  // << endmsg;
  Info << "mc file: " << mcFile ;  // << endmsg;

  std::vector<std::string> dataFileNames;
  dataFileNames.push_back(datFile);

  std::vector<std::string> mcFileNames;
  mcFileNames.push_back(mcFile);



  if(mode == "spinDensity"){

   bool cacheAmps = theAppParams->cacheAmps();
   Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
   if (cacheAmps) GlobalEnv::instance()->Channel()->Lh()->cacheAmplitudes();

   EventList eventsData;
   theAppBase.readEvents(eventsData, dataFileNames, 0, GlobalEnv::instance()->Channel()->useEvtWeight(), 0, spinDensityHist::MAX_EVENTS);

   EventList mcData;
   theAppBase.readEvents(mcData, mcFileNames, 0, false, 0, spinDensityHist::MAX_EVENTS);

   std::shared_ptr<EvtDataBaseList> pbarpEventListPtr(new EvtDataBaseList(0));
   pbarpEventListPtr->read(eventsData, mcData);

   GlobalEnv::instance()->Channel()->Lh()->setDataVec(pbarpEventListPtr->getDataVecs());
   GlobalEnv::instance()->Channel()->Lh()->setMcVec(pbarpEventListPtr->getMcVecs());

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




  EventList eventsData;
  theAppBase.readEvents(eventsData, dataFileNames, 0, GlobalEnv::instance()->Channel()->useEvtWeight());

  int maxMcEvts=eventsData.size()*ratioMcToData;

  EventList mcData;
  theAppBase.readEvents(mcData, mcFileNames, 0, false, 0, maxMcEvts-1);

  std::shared_ptr<EvtDataBaseList> pbarpEventListPtr(new EvtDataBaseList(0));

  if (mode=="plotMode"){
    std::shared_ptr<AbsHist> theHistPtr( new pbarpHist());
    theAppBase.plotMode(eventsData, mcData, pbarpEventListPtr, theHistPtr);
    return 1;
  }

  if (mode=="qaModeSimple"){
    std::shared_ptr<AbsHist> theHistPtr( new pbarpHist());
    theAppBase.qaModeSimple(eventsData, mcData, theStartparams, pbarpEventListPtr, theHistPtr, noOfFreeFitParams);
    return 1;
  }



  pbarpEventListPtr->read(eventsData, mcData);

  theLhPtr->setDataVec(pbarpEventListPtr->getDataVecs());
  theLhPtr->setMcVec(pbarpEventListPtr->getMcVecs());

  PwaFcnBase theFcn;
  Info << "\nThe parameter values are: " << "\n" << endmsg;
  GlobalEnv::instance()->fitParamsBase()->printParams(theStartparams);

  Info << "\nThe parameter errors are: " << "\n" << endmsg;
  GlobalEnv::instance()->fitParamsBase()->printParams(theErrorparams);





  if (mode=="qaMode"){
      double evtWeightSumData = pbarpEventListPtr->NoOfWeightedDataEvts();
      theAppBase.qaMode(theStartparams, evtWeightSumData, noOfFreeFitParams );
      pbarpHist theHist(theLhPtr, theStartparams);
      end= clock();
      double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
      Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;

      return 1;
  }


  if (mode=="pwa"){
    bool cacheAmps = theAppParams->cacheAmps();
    Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
    if (cacheAmps) theLhPtr->cacheAmplitudes();

    FunctionMinimum min=theAppBase.migradDefault(theFcn, upar);
    double evtWeightSumData = pbarpEventListPtr->NoOfWeightedDataEvts();
    theAppBase.printFitResult(min, theStartparams, std::cout, evtWeightSumData, noOfFreeFitParams);

    return 1;
 }

  if (mode=="evo"){

    bool cacheAmps = theAppParams->cacheAmps();
    Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
    if (cacheAmps) theLhPtr->cacheAmplitudes();

    EvoMinimizer theEvoMinimizer(theFcn, upar, GlobalEnv::instance()->parser()->evoPopulation(),
				 GlobalEnv::instance()->parser()->evoIterations());
    Info <<"start evolutionary minimizer "<< endmsg;
    std::vector<double> finalParamVec = theEvoMinimizer.Minimize();

    fitParams finalFitParams=theStartparams;
    GlobalEnv::instance()->fitParamsBase()->getFitParamVal(finalParamVec, finalFitParams);

    fitParams finalFitErrs=theErrorparams;

    std::ostringstream finalResultname;
    finalResultname << "finalResult" << outputFileNameSuffix << ".dat";

    std::ofstream theStream ( finalResultname.str().c_str() );
    GlobalEnv::instance()->fitParamsBase()->dumpParams(theStream, finalFitParams, finalFitErrs);

    return 1;
  }

 return 1;
}


