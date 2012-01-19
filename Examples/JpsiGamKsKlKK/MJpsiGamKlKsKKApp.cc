#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKParser.hh"

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKEventList.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKHist.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKReader.hh"

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKProdLh.hh"

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKData.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKFitParams.hh"


#include "PwaUtils/StreamFitParmsBase.hh"
#include "PwaUtils/PwaFcnBase.hh"
#include "PwaUtils/AbsLh.hh"

#include "Setup/PwaEnv.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "Particle/PdtParser.hh"
#include "qft++/topincludes/tensor.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"


using namespace ROOT::Minuit2;


void setErrLogMode( const JpsiGamKsKlKKParser::enErrLogMode& erlMode ) {
  switch(erlMode) {
  case JpsiGamKsKlKKParser::debug :
    ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
    break;
  case JpsiGamKsKlKKParser::trace :
    ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
    break;
  case JpsiGamKsKlKKParser::routine :
    ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
    break;
  case JpsiGamKsKlKKParser::warning :
    ErrLogger::instance()->setLevel(log4cpp::Priority::WARN);
    break;
  case JpsiGamKsKlKKParser::error :
    ErrLogger::instance()->setLevel(log4cpp::Priority::ERROR);
    break;
  case JpsiGamKsKlKKParser::alert :
    ErrLogger::instance()->setLevel(log4cpp::Priority::ALERT);
    break;
  default: 
    ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  }
}




int main(int __argc,char *__argv[]){
  clock_t start, end;
  start= clock();
  
  // Parse the command line
  static JpsiGamKsKlKKParser theAppParams(__argc, __argv);

   // Set the desired error logging mode
  setErrLogMode(theAppParams.getErrLogMode());

  std::string theCfgFile = theAppParams.getConfigFile();
  

  const std::string datFile=theAppParams.dataFile();
  const std::string mcFile=theAppParams.mcFile();
  Info << "data file: " << datFile ;  // << endmsg;
  Info << "mc file: " << mcFile ;  // << endmsg;

  ParticleTable pTable;
  PdtParser parser;
  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");
  if (!parser.parse(pdtFile, pTable)) {
    Alert << "Error: could not parse " << pdtFile ;  // << endmsg;
    exit(1);
  }

  std::vector<std::string> fileNames;

  fileNames.push_back(datFile);
  JpsiGamKsKlKKReader eventReader(fileNames, 5, 0); 
  EventList eventsData;
  eventReader.fillAll(eventsData);

  if (!eventsData.findParticleTypes(pTable))
    Warning << "could not find all particles" ;  // << endmsg;

  Info << "\nFile has " << eventsData.size() << " events. Each event has "
                  <<  eventsData.nextEvent()->size() << " final state particles.\n" ;  // << endmsg;
  eventsData.rewind();



  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = eventsData.nextEvent()) != 0 && evtCount < 20) {
    Info << "\n" 
                    << *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
                    << *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
                    << *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
                    << *(anEvent->p4(3)) << "\tm = " << anEvent->p4(3)->Mass() << "\n"
	           << *(anEvent->p4(4)) << "\tm = " << anEvent->p4(4)->Mass() << "\n"
	 ;  // << endmsg;
    ++evtCount;
  }
  eventsData.rewind();
  
  std::vector<std::string> fileNamesMc;
  fileNamesMc.push_back(mcFile);
  JpsiGamKsKlKKReader eventReaderMc(fileNamesMc, 5, 0); 
  EventList eventsMc;
  eventReaderMc.fillAll(eventsMc);
  eventsMc.rewind();

  //
  //calculate helicity angles, fill map with D-functions
  //
  boost::shared_ptr<const JpsiGamKsKlKKEventList> theJpsiGamKsKlKKEventListPtr(new JpsiGamKsKlKKEventList(eventsData, eventsMc));
  
  
  //
  //disable hypotheses, currently not in use
  //
  std::map<const std::string, bool> hypMap;
  hypMap["eta2225Hyp"]=true;
  hypMap["f02020Hyp"]=true;
  hypMap["f02020FlatteHyp"]=true;
  hypMap["f22300Hyp"]=true;
  
  
  const std::vector<std::string> disabledHyps=theAppParams.disabledHyps();
  std::vector<std::string>::const_iterator itStr;
 
  for (itStr=disabledHyps.begin(); itStr!=disabledHyps.end(); ++itStr){
    
    std::map<const std::string, bool>::const_iterator iter= hypMap.find( (*itStr) );
    if (iter !=hypMap.end()){
      hypMap[iter->first]= false;
      Info<< "hypothesis " << iter->first << " disabed" ;  // << endmsg;
    }
    else { Alert << "hypothesis " << (*itStr) << " can not be disabled";  // << endmsg;
      exit(0);
    }
  }

  boost::shared_ptr<AbsLh> theLhPtr;
  std::string startWithHyp=theAppParams.startHypo();
  
  if (startWithHyp=="production"){
    theLhPtr = boost::shared_ptr<AbsLh> (new JpsiGamKsKlKKProdLh(theJpsiGamKsKlKKEventListPtr, hypMap ));
  }
  else { 
    Alert << "start with hypothesis " << startWithHyp << " not supported!!!!" ;  // << endmsg;
    exit(1);
  }
  

  



  std::string mode=theAppParams.mode();
  if (mode=="dumpDefaultParams"){
    fitParams defaultVal;
    fitParams defaultErr;
    theLhPtr->getDefaultParams(defaultVal, defaultErr);
    std::ofstream theStreamDefault ( "defaultparams.dat");
    boost::shared_ptr<FitParamsBase> theFitParamBase=boost::shared_ptr<FitParamsBase>(new JpsiGamKsKlKKFitParams(defaultVal, defaultErr));
    theFitParamBase->dumpParams(theStreamDefault, defaultVal, defaultErr);
    return 0;      
  }



  std::string paramStreamerPath=theAppParams.fitParamFile();

  StreamFitParmsBase theParamStreamer(paramStreamerPath, boost::shared_ptr<FitParamsBase> (new JpsiGamKsKlKKFitParams()));
  fitParams theStartparams=theParamStreamer.getFitParamVal();
  fitParams theErrorparams=theParamStreamer.getFitParamErr();    

  boost::shared_ptr<FitParamsBase> theFitParamBase
    =boost::shared_ptr<FitParamsBase>(new JpsiGamKsKlKKFitParams(theStartparams, theErrorparams));
  
  if (mode=="qaMode"){
    Info << "\nThe parameter values are: " << "\n" << endmsg;
    theFitParamBase->printParams(theStartparams);

    Info << "\nThe parameter errors are: " << "\n" << endmsg;
    theFitParamBase->printParams(theErrorparams);

    double theLh=theLhPtr->calcLogLh(theStartparams);
    Info <<"theLh = "<< theLh << endmsg;

    JpsiGamKsKlKKHist theHist(theLhPtr, theStartparams);
    
    end= clock();
    double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
    Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;
    return 0;
  }
  
  
   if (mode=="pwa"){  
  PwaFcnBase theFcn(theLhPtr, theFitParamBase); 
  MnUserParameters upar;
  theFitParamBase->setMnUsrParams(upar);

  std::cout << "\n\n**************** Minuit Fit parameter **************************" << std::endl; 
  for (int i=0; i<int(upar.Params().size()); ++i){
    std::cout << upar.Name(i) << "\t" << upar.Value(i) << "\t" << upar.Error(i) << std::endl;
  }

  const std::vector<std::string> fixedParams=theAppParams.fixedParams();

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
 
  fitParams finalFitParams=theFitParamBase->getFitParamVal(finalParamVec);


  
  JpsiGamKsKlKKHist theHist(theLhPtr, finalFitParams);

  theFitParamBase->printParams(finalFitParams);
  double theLh=theLhPtr->calcLogLh(finalFitParams);
  Info <<"theLh = "<< theLh << endmsg;
  
// print and dump final fit result
   const std::vector<double> finalParamErrorVec=finalUsrParameters.Errors();
   for (size_t i=0; i<finalParamVec.size(); i++)
     {
       Info << "Value: " << finalParamVec[i] << "\t Error: " << finalParamErrorVec[i] << endmsg;
     }
   fitParams finalFitErrs=theFitParamBase->getFitParamVal(finalParamErrorVec);
   std::ofstream theStream ( "finalResult.dat");
   theFitParamBase->dumpParams(theStream, finalFitParams, finalFitErrs);
   return 0;
   }
   
  return 0;
}

