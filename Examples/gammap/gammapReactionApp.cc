//************************************************************************//
//                                                                        //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)                     //
//               Julian Pychy (julian@ep1.rub.de)                        //
//                   - Ruhr-Universität Bochum                            //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//************************************************************************//

#include <ctime>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "AppUtils/AppBase.hh"
#include "ConfigParser/gammapParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Event/EventList.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/WelcomeScreen.hh"
#include "Utils/ErrLogUtils.hh"
#include "gammapUtils/GammapChannelEnv.hh"
#include "gammapUtils/gammapBaseLh.hh"
#include "gammapUtils/gammapReaction.hh"

int main(int argc, char* argv[]){
  clock_t start = clock();

  InfoMsg << welcomeScreen << endmsg;
  InfoMsg << "Compiled " << __DATE__ << " " << __TIME__ << endmsg;

  gammapParser* parser = new gammapParser(argc, argv);
  setErrLogMode(parser->getErrLogMode());

  GlobalEnv::instance()->setup(parser);
  GlobalEnv::instance()->AddEnv(std::shared_ptr<GammapChannelEnv>(new GammapChannelEnv(parser)),
                                AbsChannelEnv::CHANNEL_GAMMAP);
  GlobalEnv::instance()->setupChannelEnvs();

  std::shared_ptr<GammapChannelEnv> channelEnv =
    std::static_pointer_cast<GammapChannelEnv>(GlobalEnv::instance()->GammapChannel());

  channelEnv->reaction()->print(std::cout);
  channelEnv->SetLh(std::shared_ptr<AbsLh>(new gammapBaseLh(channelEnv->channelID())));

  AppBase app;
  const std::string mode = parser->mode();

  if(mode == "dumpDefaultParams"){
    app.dumpDefaultParams();
    return 0;
  }

  if(mode == "dumpRandomParams"){
    app.dumpRandomParams();
    return 0;
  }

  std::shared_ptr<AbsPawianParameters> unsortedStartParams = app.streamPawianParams();
  GlobalEnv::instance()->setStartPawianParams(unsortedStartParams);
  std::shared_ptr<AbsPawianParameters> startParams = GlobalEnv::instance()->startPawianParams();

  if(mode == "gen"){
    app.generate(startParams);
    return 0;
  }

  if(parser->doScaling()){
    app.fixAllReleaseScaleParams(startParams);
  }
  else{
    app.fixParams(startParams, parser->fixedParams());
  }

  if(mode == "server" || mode == "serverGradientNum" || mode == "serverGradientNumSlow" ||
     mode == "serverAdamNum" || mode == "evoserver"){
    app.fitServerMode(startParams);
    return 0;
  }

  if(mode == "client"){
    app.fitClientMode(startParams);
    return 0;
  }

  std::shared_ptr<AbsLh> lh = channelEnv->Lh();

  std::vector<std::string> dataFileNames;
  dataFileNames.push_back(parser->dataFile());

  std::vector<std::string> mcFileNames;
  mcFileNames.push_back(parser->mcFile());

  EventList dataEvents;
  app.readEvents(dataEvents, dataFileNames, channelEnv->channelID(), channelEnv->useDataEvtWeight(),
                 0, parser->noOfDataEvts());

  EventList mcEvents;
  const int maxMcEvents = dataEvents.size() * parser->ratioMcToData();
  app.readEvents(mcEvents, mcFileNames, channelEnv->channelID(), channelEnv->useMCEvtWeight(),
                 0, maxMcEvents - 1);

  std::shared_ptr<EvtDataBaseList> eventList(new EvtDataBaseList(channelEnv->channelID()));

  if(mode == "plotMode"){
    app.plotMode(dataEvents, mcEvents, eventList);
    return 0;
  }

  if(mode == "qaModeSimple"){
    app.qaModeSimple(dataEvents, mcEvents, startParams);
    return 0;
  }

  eventList->read(dataEvents, mcEvents);
  dataEvents.removeAndDeleteEvents(0, dataEvents.size() - 1);
  mcEvents.removeAndDeleteEvents(0, mcEvents.size() - 1);

  lh->setDataVec(eventList->getDataVecs());
  lh->setMcVec(eventList->getMcVecs());

  InfoMsg << "\nThe parameter values and errors are: \n" << endmsg;
  startParams->print(std::cout);

  const double dataWeightSum = eventList->NoOfWeightedDataEvts();
  const double mcWeightSum = eventList->NoOfWeightedMcEvts();

  if(mode == "qaMode"){
    app.qaMode(startParams, dataWeightSum);
    const double cpuTime = double(clock() - start) / CLOCKS_PER_SEC;
    InfoMsg << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;
    return 0;
  }

  const bool cacheAmps = parser->cacheAmps();
  InfoMsg << "caching amplitudes enabled / disabled:\t" << cacheAmps << endmsg;
  if(cacheAmps) lh->cacheAmplitudes();

  if(mode == "pwa" || mode == "evo"){
    app.fitNonServerMode(startParams, dataWeightSum, mcWeightSum);
    return 0;
  }

  Alert << "Unsupported mode for gammapReactionApp: " << mode << endmsg;
  return 1;
}
