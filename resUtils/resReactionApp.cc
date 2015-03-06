//************************************************************************//
//									  //
//  Copyrightg 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
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
#include <iterator>
#include <memory>
#include <memory>

#include "TROOT.h"

#include "ConfigParser/resParser.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Particle/PdtParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/AbsLh.hh"
#include "FitParams/FitParColBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "FitParams/StreamFitParColBase.hh"
#include "MinFunctions/PwaFcnBase.hh"
#include "PwaUtils/WelcomeScreen.hh"
#include "FitParams/PwaCovMatrix.hh"
#include "PwaUtils/WaveContribution.hh"
#include "AppUtils/AppBase.hh"
#include "PwaUtils/NetworkClient.hh"
#include "PwaUtils/NetworkServer.hh"
#include "MinFunctions/PwaFcnServer.hh"
#include "MinFunctions/EvoMinimizer.hh"

#include "Utils/PawianCollectionUtils.hh"
#include "Utils/ErrLogUtils.hh"
#include "resUtils/ResChannelEnv.hh"
#include "resUtils/resReaction.hh"
#include "resUtils/resBaseLh.hh"

#include "Event/EventReaderDefault.hh"
#include "resUtils/resHist.hh"

#include "Event/Event.hh"
#include "Event/EventList.hh"

#include "FitParams/AbsPawianParameters.hh"
#include "FitParams/MnPawianParameters.hh"
#include "FitParams/ParamFactory.hh"

int main(int __argc,char *__argv[]){
  clock_t start, end;
  start= clock();

  Info << welcomeScreen << endmsg;
  Info << "Compiled " << __DATE__ << " " << __TIME__ << endmsg;

  // Disable output buffering
  setvbuf(stdout, NULL, _IONBF, 0);

  // Parse the command line
  resParser* theAppParams=new resParser(__argc, __argv);

  // Set the desired error logging mode
  setErrLogMode(theAppParams->getErrLogMode());

  // Setup the global environment and add the primary res channel
  GlobalEnv::instance()->setup(theAppParams);
  GlobalEnv::instance()->AddEnv(std::shared_ptr<ResChannelEnv>(new ResChannelEnv(theAppParams)), AbsChannelEnv::CHANNEL_RES);

  // Print the primary reaction
  std::shared_ptr<resReaction> theResReaction=std::static_pointer_cast<ResChannelEnv>(GlobalEnv::instance()->ResChannel())->reaction();
  theResReaction->print(std::cout);

  //
  //! resReactionApp does not support coupled channel analyses yet
  //

 std::string mode=theAppParams->mode();


 std::string prodFormalism=theAppParams->productionFormalism();
 std::shared_ptr<AbsLh> theLhPtr;
 theLhPtr = std::shared_ptr<AbsLh>(new resBaseLh(0));

 AppBase theAppBase;

 if (mode=="dumpDefaultParams"){
     theAppBase.dumpDefaultParams();
      return 0;
 }


 std::string paramStreamerPath=theAppParams->fitParamFile();
  std::string outputFileNameSuffix= GlobalEnv::instance()->outputFileNameSuffix();
  StreamFitParmsBase theParamStreamer(paramStreamerPath);
  fitParCol theStartparams=theParamStreamer.getFitParamVal();
  fitParCol theErrorparams=theParamStreamer.getFitParamErr();

  if (mode=="gen"){
    theAppBase.generate(theStartparams);
    return 0;
  }


  const std::string datFile=theAppParams->dataFile();
  const std::string mcFile=theAppParams->mcFile();
  Info << "data file: " << datFile ;  // << endmsg;
  Info << "mc file: " << mcFile ;  // << endmsg;

  std::vector<std::string> dataFileNames;
  dataFileNames.push_back(datFile);

  std::vector<std::string> mcFileNames;
  mcFileNames.push_back(mcFile);

  bool withEvtWeight=theAppParams->useEvtWeight();
  Info << "EvtWeight: " << withEvtWeight << endmsg;

  std::shared_ptr<AbsPawianParameters> upar=ParamFactory::instance()->getParametersPointer("Minuit2");
  GlobalEnv::instance()->fitParColBase()->setAbsPawianParams(upar, theStartparams, theErrorparams);

  std::cout << "\n\n**************** Fit parameter **************************" << std::endl;
  for (int i=0; i<int(upar->Params().size()); ++i){
    std::cout << upar->Name(i) << "\t" << upar->Value(i) << "\t" << upar->Error(i) << std::endl;
  }

  const std::vector<std::string> fixedParams=theAppParams->fixedParams();
  const unsigned int noOfFreeFitParams = upar->Params().size()-fixedParams.size();

  if(mode == "client"){

  bool cacheAmps = theAppParams->cacheAmps();
  Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
  if (cacheAmps) theLhPtr->cacheAmplitudes();

  std::ostringstream portStringStream;
  portStringStream << theAppParams->serverPort();

  NetworkClient theClient(theAppParams->serverAddress(), portStringStream.str());
  if(!theClient.Login())
    return 0;


  EventList eventsDataClient;
  theAppBase.readEvents(eventsDataClient, dataFileNames, withEvtWeight, theClient.GetEventLimits()[0], theClient.GetEventLimits()[1]);

  EventList mcDataClient;
  theAppBase.readEvents(mcDataClient, mcFileNames, withEvtWeight, theClient.GetEventLimits()[2], theClient.GetEventLimits()[3]);

  std::shared_ptr<EvtDataBaseList> epemEventListPtr(new EvtDataBaseList(0));
  epemEventListPtr->read(eventsDataClient, mcDataClient);

  theLhPtr->setDataVec(epemEventListPtr->getDataVecs());
  theLhPtr->setMcVec(epemEventListPtr->getMcVecs());

  theAppBase.calcAndSendClientLh(theClient, theStartparams, 0);

  return 1;
 }

  EventList eventsData;
  theAppBase.readEvents(eventsData, dataFileNames, withEvtWeight);

  int ratioMcToData=theAppParams->ratioMcToData();
  int maxMcEvts=eventsData.size()*ratioMcToData;
  EventList mcData;
  theAppBase.readEvents(mcData, mcFileNames, withEvtWeight, 0, maxMcEvts-1);

  std::shared_ptr<EvtDataBaseList> eventListPtr(new EvtDataBaseList(0));
  eventListPtr->read(eventsData, mcData);

  theLhPtr->setDataVec(eventListPtr->getDataVecs());
  theLhPtr->setMcVec(eventListPtr->getMcVecs());

  if(mode == "server"){
    theAppBase.fixParams(upar,fixedParams);

    std::map<short, std::tuple<long, double, long> > numEventMap;
    numEventMap[0] = std::tuple<long, double,long>(eventsData.size(), eventsData.size(), mcData.size());

    std::shared_ptr<NetworkServer> theServer(new NetworkServer(theAppParams->serverPort(),
 							       theAppParams->noOfClients(),
 							       numEventMap, theAppParams->clientNumberWeights()));


    PwaFcnServer theFcnServer(theServer);
    theServer->WaitForFirstClientLogin();

    FunctionMinimum min=theAppBase.migradDefault(theFcnServer, upar);

    theServer->BroadcastClosingMessage();
    Info << "Closing server." << endmsg;

    double evtWeightSumData = eventListPtr->NoOfWeightedDataEvts();
    theAppBase.printFitResult(min, theStartparams, std::cout, evtWeightSumData, noOfFreeFitParams);

    return 1;
 }

 if(mode == "evoserver"){
   theAppBase.fixParams(upar,fixedParams);

   std::map<short, std::tuple<long, double, long> > numEventMap;
   numEventMap[0] = std::tuple<long, double,long>(eventsData.size(), eventsData.size(), mcData.size());

    std::shared_ptr<NetworkServer> theServer(new NetworkServer(theAppParams->serverPort(),
 							       theAppParams->noOfClients(),
 							       numEventMap, theAppParams->clientNumberWeights()));

   PwaFcnServer theFcnServer(theServer);
   theServer->WaitForFirstClientLogin();

   EvoMinimizer theEvoMinimizer(theFcnServer, upar, GlobalEnv::instance()->parser()->evoPopulation(),
                                GlobalEnv::instance()->parser()->evoIterations());
   Info <<"start evolutionary minimizer "<< endmsg;
   std::vector<double> finalParamVec = theEvoMinimizer.Minimize();

   theServer->BroadcastClosingMessage();
   Info << "Closing server." << endmsg;

   fitParCol finalFitParams=theStartparams;
   GlobalEnv::instance()->fitParColBase()->getFitParamVal(finalParamVec, finalFitParams);

   fitParCol finalFitErrs=theErrorparams;

   std::ostringstream finalResultname;
   finalResultname << "finalResult" << outputFileNameSuffix << ".dat";

   std::ofstream theStream ( finalResultname.str().c_str() );
   GlobalEnv::instance()->fitParColBase()->dumpParams(theStream, finalFitParams, finalFitErrs);
   return 1;
 }


  PwaFcnBase theFcn;

  if (mode=="qaMode"){
    double evtWeightSumData = eventListPtr->NoOfWeightedDataEvts();
    theAppBase.qaMode(theStartparams, evtWeightSumData, noOfFreeFitParams );
    resHist theHist;
    theHist.fillFromLhData(theLhPtr, theStartparams);
    end= clock();
    double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
    Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;

    return 1;
  }


  if (mode=="pwa"){
    bool cacheAmps = theAppParams->cacheAmps();
    Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
    if (cacheAmps) theLhPtr->cacheAmplitudes();

    theAppBase.fixParams(upar, fixedParams);

    FunctionMinimum min=theAppBase.migradDefault(theFcn, upar);

    double evtWeightSumData = eventListPtr->NoOfWeightedDataEvts();
    theAppBase.printFitResult(min, theStartparams, std::cout, evtWeightSumData, noOfFreeFitParams);

    return 1;
 }

}


