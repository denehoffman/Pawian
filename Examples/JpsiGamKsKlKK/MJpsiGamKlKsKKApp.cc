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

#include "Examples/JpsiGamKsKlKK/AbsJpsiGamKsKlKKLh.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKProdLh.hh"
//#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKBaseLh.hh"


#include "Examples/JpsiGamKsKlKK/MJpsiGamKsKlKKFcn.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKData.hh"

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKStreamFitParams.hh"

#include "Setup/PwaEnv.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
// #include "Event/CBElsaReader.hh"
#include "Particle/PdtParser.hh"
#include "qft++/topincludes/tensor.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"


//#include "Minuit2/MnUserTransformation.h"
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

  // Parse the command line
  static JpsiGamKsKlKKParser theAppParams(__argc, __argv);

   // Set the desired error logging mode
  setErrLogMode(theAppParams.getErrLogMode());

  std::string theCfgFile = theAppParams.getConfigFile();
  Info << "The path to config file is " << theCfgFile << "\n" ;  // << endmsg;    

  std::string paramStreamerPath=theAppParams.fitParamFile();

  JpsiGamKsKlKKStreamFitParams theParamStreamer(paramStreamerPath);
  paramGamKsKlKK theStartparams=theParamStreamer.getFitParamVal();
  paramGamKsKlKK theErrorparams=theParamStreamer.getFitParamErr();

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
  hypMap["K0_1430HypBase"]=true;
  hypMap["K1_1410HypBase"]=true;
  hypMap["K2_1430HypBase"]=true;
  hypMap["K1_1680_HypBase"]=true;
  hypMap["KKPi_HypBase"]=true;

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

  boost::shared_ptr<AbsJpsiGamKsKlKKLh> theJpsiGamKsKlKKLhPtr;
  std::string startWithHyp=theAppParams.startHypo();
  
  if (startWithHyp=="production"){
    theJpsiGamKsKlKKLhPtr = boost::shared_ptr<AbsJpsiGamKsKlKKLh> (new JpsiGamKsKlKKProdLh(theJpsiGamKsKlKKEventListPtr) );
  }
//   else if (startWithHyp=="base") 
//     theJpsiGamKsKlKKLhPtr = boost::shared_ptr<AbsJpsiGamKsKlKKLh> (new JpsiGamKsKlKKBaseLh(theJpsiGamKsKlKKEventListPtr, hypMap));
  
  else { 
    Alert << "start with hypothesis " << startWithHyp << " not supported!!!!" ;  // << endmsg;
    exit(1);
  }
  

  bool qaMode=theAppParams.qaMode();
  std::cout << "qaMode: " << qaMode << std::endl;

  if (qaMode){
    theJpsiGamKsKlKKLhPtr->printCurrentFitResult(theStartparams);
    double theLh=theJpsiGamKsKlKKLhPtr->calcLogLh(theStartparams);
    Info <<"theLh = "<< theLh ;  // << endmsg;

    JpsiGamKsKlKKHist theJpsiGamKsKlKKHist(theJpsiGamKsKlKKLhPtr, theStartparams);
    return 0;
  }
  

   MJpsiGamKsKlKKFcn mJpsiGamKsKlKKFcn(theJpsiGamKsKlKKLhPtr);

   MnUserParameters upar;
   theJpsiGamKsKlKKLhPtr->setMnUsrParams(upar, theStartparams, theErrorparams);

  const std::vector<std::string> fixedParams=theAppParams.fixedParams();

  std::vector<std::string>::const_iterator itFix;
  for (itFix=fixedParams.begin(); itFix!=fixedParams.end(); ++itFix){
    upar.Fix( (*itFix) );
  }

  
  MnMigrad migrad(mJpsiGamKsKlKKFcn, upar);
  vector<ROOT::Minuit2::MinuitParameter> vec = migrad.MinuitParameters();
  vector<ROOT::Minuit2::MinuitParameter>::iterator iter;
  for(iter=vec.begin(); iter!=vec.end();++iter){
    Info << "Parameter " << (*iter).Name() << " " << (*iter).Value() ;  // << endmsg;
  } 
  theJpsiGamKsKlKKLhPtr->printCurrentFitResult(theStartparams);
  


  Info <<"start migrad ";  // << endmsg;
  FunctionMinimum min = migrad();

  if(!min.IsValid()) {
    //try with higher strategy
    Info <<"FM is invalid, try with strategy = 2.";  // << endmsg;
    MnMigrad migrad2(mJpsiGamKsKlKKFcn, min.UserState(), MnStrategy(2));
    min = migrad2();
  }

  MnUserParameters finalUsrParameters=min.UserParameters();
  const std::vector<double> finalParamVec=finalUsrParameters.Params();
 
  paramGamKsKlKK finalFitParams;
  theJpsiGamKsKlKKLhPtr->setFitParamVal(finalFitParams, finalParamVec);


  
  JpsiGamKsKlKKHist theJpsiGamKsKlKKHist(theJpsiGamKsKlKKLhPtr, finalFitParams);
  
  Info << "Final fit results: ------------------------------------------";
  theJpsiGamKsKlKKLhPtr->printCurrentFitResult(finalFitParams);
  double theLh=theJpsiGamKsKlKKLhPtr->calcLogLh(finalFitParams);
  Info <<"theLh = "<< theLh ;  // << endmsg;
  
  // print final fit result
  
  const std::vector<double> finalParamErrorVec=finalUsrParameters.Errors();
  for (size_t i=0; i<finalParamVec.size(); i++)
    {
      Info << "Value: " << finalParamVec[i] << "\t Error: " << finalParamErrorVec[i] ;  // << endmsg;
    }



    return 0;
}

