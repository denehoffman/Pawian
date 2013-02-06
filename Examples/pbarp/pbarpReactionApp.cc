#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <map>
#include <iterator>
#include <boost/shared_ptr.hpp>

#include "TROOT.h"

#include "pbarpUtils/pbarpParser.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Particle/PdtParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include "pbarpUtils/pbarpStatesLS.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/StreamFitParmsBase.hh"
#include "PwaUtils/PwaFcnBase.hh"
#include "PwaUtils/PwaCovMatrix.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/ErrLogUtils.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "pbarpUtils/pbarpReaction.hh"
#include "pbarpUtils/pbarpBaseLh.hh"
#include "pbarpUtils/pbarpHeliLh.hh"
#include "pbarpUtils/pbarpCanoLh.hh"
#include "pbarpUtils/pbarpEvtReader.hh"
#include "PwaUtils/EvtDataBaseList.hh"
//#include "pbarpUtils/pbarpEventList.hh"
#include "pbarpUtils/pbarpHist.hh"
#include "pbarpUtils/spinDensityHist.hh"
#include "Event/Event.hh"
#include "Event/EventList.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnScan.h"


#ifdef _OPENMP
#include <omp.h>
#endif


int main(int __argc,char *__argv[]){
  clock_t start, end;
  start= clock();

  // Disable output buffering
  setvbuf(stdout, NULL, _IONBF, 0);

  // Parse the command line
  pbarpParser* theAppParams=new pbarpParser(__argc, __argv);

  // Set the desired error logging mode
  setErrLogMode(theAppParams->getErrLogMode());
  
#ifdef _OPENMP
  const int noOfThreads=theAppParams->noOfThreads();
  omp_set_num_threads(noOfThreads);
#endif


  pbarpEnv::instance()->setup(theAppParams);

  boost::shared_ptr<pbarpReaction> thepbarpReaction=pbarpEnv::instance()->reaction();

  thepbarpReaction->print(std::cout);

  // boost::shared_ptr<pbarpDataBaseList> thepbarbDataBaseListPtr(new pbarpDataBaseList()); 
  // boost::shared_ptr<AbsLh> theLhPtr(new pbarpBaseLh(thepbarbDataBaseListPtr));

  std::string mode=theAppParams->mode();

  boost::shared_ptr<FitParamsBase> theFitParamBase=boost::shared_ptr<FitParamsBase>(new FitParamsBase());



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


  //  boost::shared_ptr<pbarpEventList> pbarpEventListPtr(new pbarpEventList());
  boost::shared_ptr<EvtDataBaseList> pbarpEventListPtr(new EvtDataBaseList(pbarpEnv::instance()));
  pbarpEventListPtr->ratioMcToData(theAppParams->ratioMcToData());
  pbarpEventListPtr->read(eventsData, mcData);
  // pbarpEventListPtr->read4Vecs();

  std::string prodFormalism=theAppParams->productionFormalism();
  boost::shared_ptr<AbsLh> theLhPtr;
  if(prodFormalism=="Cano") theLhPtr=boost::shared_ptr<AbsLh>(new pbarpCanoLh(pbarpEventListPtr));
  else if(prodFormalism=="Heli") theLhPtr=boost::shared_ptr<AbsLh>(new pbarpHeliLh(pbarpEventListPtr));
  else {
    Alert << "prodFormalism\t" << prodFormalism << "\tdoesn't exist!!!" << endmsg;
    exit(1);
  }

  if (mode=="dumpDefaultParams"){
    fitParams defaultVal;
    fitParams defaultErr;
    theLhPtr->getDefaultParams(defaultVal, defaultErr);

    std::stringstream defaultparamsname;
    defaultparamsname << "defaultparams" << pbarpEnv::instance()->outputFileNameSuffix() << ".dat";
    std::ofstream theStreamDefault ( defaultparamsname.str().c_str() );
    
    theFitParamBase->dumpParams(theStreamDefault, defaultVal, defaultErr);
    return 0;
  }


  std::string paramStreamerPath=theAppParams->fitParamFile();
  std::string outputFileNameSuffix= pbarpEnv::instance()->outputFileNameSuffix();
  StreamFitParmsBase theParamStreamer(paramStreamerPath, theLhPtr);
  fitParams theStartparams=theParamStreamer.getFitParamVal();
  fitParams theErrorparams=theParamStreamer.getFitParamErr();

  PwaFcnBase theFcn(theLhPtr, theFitParamBase, outputFileNameSuffix);
  MnUserParameters upar;
  theFitParamBase->setMnUsrParams(upar, theStartparams, theErrorparams);
  
  std::cout << "\n\n**************** Minuit Fit parameter **************************" << std::endl;
  for (int i=0; i<int(upar.Params().size()); ++i){
    std::cout << upar.Name(i) << "\t" << upar.Value(i) << "\t" << upar.Error(i) << std::endl;
  }

  const std::vector<std::string> fixedParams=theAppParams->fixedParams();  
  const unsigned int noOfFreeFitParams = upar.Params().size()-fixedParams.size();

  if (mode=="qaMode"){
    Info << "\nThe parameter values are: " << "\n" << endmsg;
    theFitParamBase->printParams(theStartparams);
    
    Info << "\nThe parameter errors are: " << "\n" << endmsg;
    theFitParamBase->printParams(theErrorparams);
    
    double theLh=theLhPtr->calcLogLh(theStartparams);
    Info <<"theLh = "<< theLh << endmsg;

    pbarpHist theHist(theLhPtr, theStartparams);

    double evtWeightSumData = pbarpEventListPtr->NoOfWeightedDataEvts();
    double BICcriterion=2.*theLh+noOfFreeFitParams*log(evtWeightSumData);
    double AICcriterion=2.*theLh+2.*noOfFreeFitParams;
    double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)/(evtWeightSumData-noOfFreeFitParams-1);
    
    Info << "noOfFreeFitParams:\t" <<noOfFreeFitParams;
    Info << "evtWeightSumData:\t" <<evtWeightSumData; 
    Info << "BIC:\t" << BICcriterion << endmsg;
    Info << "AIC:\t" << AICcriterion << endmsg;
    Info << "AICc:\t" << AICccriterion << endmsg;
    
    std::ostringstream qaSummaryFileName;
    qaSummaryFileName << "qaSummary" << outputFileNameSuffix << ".dat";

    std::ofstream theQaStream ( qaSummaryFileName.str().c_str() );
    theQaStream << "BIC\t" << BICcriterion << "\n";
    theQaStream << "AICa\t" << AICcriterion << "\n";
    theQaStream << "AICc\t" << AICccriterion << "\n";
    theQaStream << "logLh\t" << theLh << "\n";
    theQaStream << "free parameter\t" << noOfFreeFitParams << "\n";
    theQaStream.close();
    
    end= clock();
    double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
    Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;
    
    return 1;    
  }
  

  if (mode=="pwa"){
    bool cacheAmps = theAppParams->cacheAmps();
    Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
    if (cacheAmps) theLhPtr->cacheAmplitudes();
    std::vector<std::string>::const_iterator itFix;
    for (itFix=fixedParams.begin(); itFix!=fixedParams.end(); ++itFix){
      upar.Fix( (*itFix) );
    }

    MnMigrad migrad(theFcn, upar);
    Info <<"start migrad "<< endmsg;
    FunctionMinimum min = migrad();
    
    if(!min.IsValid()) {
      //try with higher strategy
      Info <<"FM is invalid, try with strategy = 2."<< endmsg;
      MnMigrad migrad2(theFcn, min.UserState(), MnStrategy(2));
      min = migrad2();
    }
    
    MnUserParameters finalUsrParameters=min.UserParameters();
    const std::vector<double> finalParamVec=finalUsrParameters.Params();
    fitParams finalFitParams=theStartparams;
    theFitParamBase->getFitParamVal(finalParamVec, finalFitParams);

    theFitParamBase->printParams(finalFitParams);
    double theLh=theLhPtr->calcLogLh(finalFitParams);
    Info <<"theLh = "<< theLh << endmsg;
    
    
    const std::vector<double> finalParamErrorVec=finalUsrParameters.Errors();
    fitParams finalFitErrs=theErrorparams;
    theFitParamBase->getFitParamVal(finalParamErrorVec, finalFitErrs);
    
    std::ostringstream finalResultname;
    finalResultname << "finalResult" << outputFileNameSuffix << ".dat";

    std::ofstream theStream ( finalResultname.str().c_str() );
    theFitParamBase->dumpParams(theStream, finalFitParams, finalFitErrs);
    
    MnUserCovariance theCovMatrix = min.UserCovariance();
    std::cout  << min << std::endl;

    std::ostringstream serializationFileName;
    serializationFileName << "serializedOutput" << outputFileNameSuffix << ".dat";
    std::ofstream serializationStream(serializationFileName.str().c_str());
    boost::archive::text_oarchive boostOutputArchive(serializationStream);

    if(min.IsValid()){
       PwaCovMatrix thePwaCovMatrix(theCovMatrix, finalUsrParameters, finalFitParams);
       boostOutputArchive << thePwaCovMatrix;
    }

    return 0;
 }

 if(mode == "spinDensity"){

    std::string serializationFileName = pbarpEnv::instance()->serializationFileName();
    std::ifstream serializationStream(serializationFileName.c_str());

    if(!serializationStream.is_open()){
       Alert << "Could not open serialization file." << endmsg;
       return 0;
    }

    boost::archive::text_iarchive boostInputArchive(serializationStream);

    PwaCovMatrix thePwaCovMatrix;
    boostInputArchive >> thePwaCovMatrix;

    spinDensityHist theSpinDensityHists(theLhPtr, theStartparams, thePwaCovMatrix);
 }
}     
 
  
