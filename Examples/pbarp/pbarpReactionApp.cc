#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <map>
#include <iterator>

#include "Examples/pbarp/pbarpParser.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Particle/PdtParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/pbarpStatesLS.hh"
#include "PwaUtils/AbsLhNew.hh"
#include "PwaUtils/FitParamsBaseNew.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/ErrLogUtils.hh"
#include "Examples/pbarp/IsobarDecay.hh"
#include "Examples/pbarp/IsobarDecayList.hh"
#include "Examples/pbarp/pbarpEnv.hh"
#include "Examples/pbarp/pbarpReaction.hh"
#include "Examples/pbarp/pbarpDataBaseList.hh"
#include "Examples/pbarp/pbarpBaseLh.hh"
#include "Examples/pbarp/pbarpEvtReader.hh"
#include "Examples/pbarp/pbarpEventList.hh"
#include "Event/Event.hh"
#include "Event/EventList.hh"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int __argc,char *__argv[]){
  clock_t start, end;
  start= clock();

  // Parse the command line
  static pbarpParser theAppParams(__argc, __argv);

  // Set the desired error logging mode
  setErrLogMode(theAppParams.getErrLogMode());
  
#ifdef _OPENMP
  const int noOfThreads=theAppParams.noOfThreads();
  omp_set_num_threads(noOfThreads);
#endif


  pbarpEnv::instance()->setup(theAppParams);

  boost::shared_ptr<pbarpReaction> thepbarpReaction=pbarpEnv::instance()->reaction();

  thepbarpReaction->print(std::cout);

  boost::shared_ptr<pbarpDataBaseList> thepbarbDataBaseListPtr(new pbarpDataBaseList()); 
  boost::shared_ptr<AbsLhNew> theLhPtr(new pbarpBaseLh(thepbarbDataBaseListPtr));

  std::string mode=theAppParams.mode();

  boost::shared_ptr<FitParamsBaseNew> theFitParamBase=boost::shared_ptr<FitParamsBaseNew>(new FitParamsBaseNew());

  if (mode=="dumpDefaultParams"){
    fitParamsNew defaultVal;
    fitParamsNew defaultErr;
    theLhPtr->getDefaultParams(defaultVal, defaultErr);
    std::string defaultparamsname = "defaultparams.dat";
    std::ofstream theStreamDefault ( defaultparamsname.c_str() );
    
    theFitParamBase->dumpParams(theStreamDefault, defaultVal, defaultErr);
    return 0;
  }

  const std::string datFile=theAppParams.dataFile();
  const std::string mcFile=theAppParams.mcFile();
  Info << "data file: " << datFile ;  // << endmsg;
  Info << "mc file: " << mcFile ;  // << endmsg;
  
  std::vector<std::string> dataFileNames;
  dataFileNames.push_back(datFile);

  std::vector<std::string> mcFileNames;
  mcFileNames.push_back(mcFile);  

  bool withEvtWeight=theAppParams.useEvtWeight();
  Info << "EvtWeight: " << withEvtWeight << endmsg;  

  
  int noFinalStateParticles=pbarpEnv::instance()->noFinalStateParticles();  
  pbarpEvtReader eventReaderData(dataFileNames, noFinalStateParticles, 0, withEvtWeight);

  EventList eventsData;
  eventReaderData.fillAll(eventsData);
  
  Info  << "\nFile has " << eventsData.size() << " events. Each event has "
        <<  eventsData.nextEvent()->size() << " final state particles.\n" ;  // << endmsg;
  eventsData.rewind();

  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = eventsData.nextEvent()) != 0 && evtCount < 10) {
    Info        << "\n";
    for(int i=0; i<noFinalStateParticles; ++i){
      Info        << (*anEvent->p4(i)) << "\tm = " << anEvent->p4(i)->Mass() << "\n";
    }
    Info        << "\n" << endmsg;
    ;  // << endmsg;
    ++evtCount;
  }
  eventsData.rewind();

  pbarpEvtReader eventReaderMc(mcFileNames, noFinalStateParticles, 0, false);

  EventList mcData;
  eventReaderMc.fillAll(mcData);
 Info  << "\nFile has " << mcData.size() << " events. Each event has "
        <<  mcData.nextEvent()->size() << " final state particles.\n" ;  // << endmsg;
  mcData.rewind();

  evtCount = 0;
  while ((anEvent = mcData.nextEvent()) != 0 && evtCount < 10) {
    Info        << "\n";
    for(int i=0; i<noFinalStateParticles; ++i){
      Info        << (*anEvent->p4(i)) << "\tm = " << anEvent->p4(i)->Mass() << "\n";
    }
    Info        << "\n" << endmsg;
    ;  // << endmsg;
    ++evtCount;
  }
  mcData.rewind();


  boost::shared_ptr<pbarpEventList> pbarpEventListPtr(new pbarpEventList());
  pbarpEventListPtr->ratioMcToData(theAppParams.ratioMcToData());
  pbarpEventListPtr->read(eventsData, mcData);
  // pbarpEventListPtr->read4Vecs();
}     
 
  


