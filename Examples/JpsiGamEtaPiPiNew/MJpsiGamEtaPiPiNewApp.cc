#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <boost/shared_ptr.hpp>

#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiParser.hh"
#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiReader.hh"
#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiHistNew.hh"
#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiEventListNew.hh"
//#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiStates.hh"
#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiProdLhNew.hh"

#include "PwaUtils/EvtDataBaseListNew.hh"
#include "PwaUtils/FitParamsBaseNew.hh"
#include "PwaUtils/StreamFitParmsBaseNew.hh"
#include "PwaUtils/PwaFcnBaseNew.hh"
#include "PwaUtils/AbsLhNew.hh"

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
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnScan.h"

using namespace ROOT::Minuit2;


void setErrLogMode( const JpsiGamEtaPiPiParser::enErrLogMode& erlMode ) {
  std::cout << erlMode << std::endl;
	switch(erlMode) {
	case JpsiGamEtaPiPiParser::debug :
		ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
		break;
	case JpsiGamEtaPiPiParser::trace :
		ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
		break;
	case JpsiGamEtaPiPiParser::routine :
		ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
		break;
	case JpsiGamEtaPiPiParser::warning :
		ErrLogger::instance()->setLevel(log4cpp::Priority::WARN);
		break;
	case JpsiGamEtaPiPiParser::error :
		ErrLogger::instance()->setLevel(log4cpp::Priority::ERROR);
		break;
	case JpsiGamEtaPiPiParser::alert :
		ErrLogger::instance()->setLevel(log4cpp::Priority::ALERT);
		break;
	default:
		ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
	}
}

int main(int __argc,char *__argv[]){
  setvbuf(stdout, NULL, _IONBF, 0); // Outputbuffer
  clock_t start, end;
  start= clock();
  
  // Parse the command line
  static JpsiGamEtaPiPiParser theAppParams(__argc, __argv);
  
  // Set the desired error logging mode
  setErrLogMode(theAppParams.getErrLogMode());
  
#ifdef _OPENMP
  const int noOfThreads=theAppParams.noOfThreads();
  omp_set_num_threads(noOfThreads);
#endif
  
  
  std::string theCfgFile = theAppParams.getConfigFile();
  std::string jobOption = theAppParams.getjobOption();
  
  const std::string datFile=theAppParams.dataFile();
  const std::string mcFile=theAppParams.mcFile();
  const std::string sumFile=theAppParams.sumFile();
  Info << "data file: " << datFile ;  // << endmsg;
  Info << "mc file: " << mcFile ;  // << endmsg;
  Info << "sum file: " << sumFile ;  // << endmsg;
  
  ParticleTable pTable;
  PdtParser parser;
  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");
  if (!parser.parse(pdtFile, pTable)) {
    Alert << "Error: could not parse " << pdtFile ;  // << endmsg;
    exit(1);
  }
  
  std::pair<double, double> massRange = theAppParams.massRange();
  Info  << "Mass range: " << massRange.first << " " << massRange.second ;
  
  std::vector<std::string> fileNames;
  fileNames.push_back(datFile);
  JpsiGamEtaPiPiReader eventReader(fileNames, 4, 0);
  EventList eventsData;
  eventReader.fillMassRange(eventsData, massRange );
  
  if (!eventsData.findParticleTypes(pTable))
    Warning << "could not find all particles" ;  // << endmsg;
  
  Info 	<< "\nFile has " << eventsData.size() << " events. Each event has "
	<<  eventsData.nextEvent()->size() << " final state particles.\n" ;  // << endmsg;
  eventsData.rewind();
  
  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = eventsData.nextEvent()) != 0 && evtCount < 1) {
    Info	<< "\n"
		<< *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
		<< *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
		<< *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
		<< *(anEvent->p4(3)) << "\tm = " << anEvent->p4(3)->Mass() << "\n"
      ;  // << endmsg;
    ++evtCount;
  }
  eventsData.rewind();
  
  
  std::vector<std::string> fileNamesMc;
  fileNamesMc.push_back(mcFile);
  JpsiGamEtaPiPiReader eventReaderMc(fileNamesMc, 4, 0);
  EventList eventsMc;
  eventReaderMc.fillMassRange(eventsMc, massRange);
  eventsMc.rewind();
  
  //
  //calculate helicity angles, fill map with D-functions
  //
  
  boost::shared_ptr<const JpsiGamEtaPiPiEventListNew> theJpsiGamEtaPiPiEventListPtr(new JpsiGamEtaPiPiEventListNew(eventsData, eventsMc));
  
  std::string mode=theAppParams.mode();
  std::cout << "Mode: " << mode << std::endl;
  if (mode=="plotmode"){
    JpsiGamEtaPiPiHistNew theHist(theJpsiGamEtaPiPiEventListPtr,theAppParams.massRange(), jobOption.c_str());
    return 0;
  }
  
  //
  //retrieve  hypotheses
  //
  
  boost::shared_ptr<JpsiGamEtaPiPiStates> jpsiGamEtaPiPiStatesPtr(new JpsiGamEtaPiPiStates());
  const std::vector<std::string> hypVec=theAppParams.enabledHyps();
  boost::shared_ptr<AbsLhNew> theLhPtr;
  
  std::string startWithHyp=theAppParams.startHypo();
  
  if (startWithHyp=="production"){
    theLhPtr = boost::shared_ptr<AbsLhNew> (new JpsiGamEtaPiPiProdLhNew(theJpsiGamEtaPiPiEventListPtr, hypVec, jpsiGamEtaPiPiStatesPtr));
  }
  else {
    Alert << "start with hypothesis " << startWithHyp << " not supported!!!!" ;  // << endmsg;
    exit(1);
  }
  
  boost::shared_ptr<FitParamsBaseNew> theFitParamBase=boost::shared_ptr<FitParamsBaseNew>(new FitParamsBaseNew());
  
  if (mode=="dumpDefaultParams"){
    fitParamsNew defaultVal;
    fitParamsNew defaultErr;
    theLhPtr->getDefaultParams(defaultVal, defaultErr);
    std::string defaultparamsname = "defaultparams" + jobOption + ".dat";
    std::ofstream theStreamDefault ( defaultparamsname.c_str() );
    //	  std::ofstream theStreamDefault ( "defaultparams.dat");
    
    theFitParamBase->dumpParams(theStreamDefault, defaultVal, defaultErr);
    return 0;
  }
  
  std::string paramStreamerPath=theAppParams.fitParamFile();
  StreamFitParmsBaseNew theParamStreamer(paramStreamerPath, theLhPtr);
  fitParamsNew theStartparams=theParamStreamer.getFitParamVal();
  fitParamsNew theErrorparams=theParamStreamer.getFitParamErr();
  
  if (mode=="qaMode"){
    Info << "\nThe parameter values are: " << "\n" << endmsg;
    theFitParamBase->printParams(theStartparams);
    Info << "\nThe parameter errors are: " << "\n" << endmsg;
    theFitParamBase->printParams(theErrorparams);
    
    double theLh=theLhPtr->calcLogLh(theStartparams);
    Info <<"theLh = "<< theLh << endmsg;
    
    JpsiGamEtaPiPiHistNew theHist(theLhPtr, theStartparams,theAppParams.massRange(), jobOption.c_str());
    theHist.PrintToPDF(jobOption);
    theHist.SaveToROOT();

    // 		if(theAppParams.massIndependentFit()){
    // 			//calculate intensity contributions
    // 			//std::ofstream theStream ( "componentIntensity.dat");
    // 			//theProdLh->dumpComponentIntensity( theStream, theStartparams, theErrorMatrix );
    // 		}
    
    end= clock();
    double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
    Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;
    return 0;
  }

  if (mode=="pwa"){
    PwaFcnBaseNew theFcn(theLhPtr, theFitParamBase, jobOption);
    MnUserParameters upar;
    theFitParamBase->setMnUsrParams(upar, theStartparams, theErrorparams);
    
    std::cout << "\n\n**************** Minuit Fit parameter **************************" << std::endl;
    for (int i=0; i<int(upar.Params().size()); ++i){
      std::cout << upar.Name(i) << "\t" << upar.Value(i) << "\t" << upar.Error(i) << std::endl;
    }
    
    const std::vector<std::string> fixedParams=theAppParams.fixedParams();
    
    std::vector<std::string>::const_iterator itFix;
    for (itFix=fixedParams.begin(); itFix!=fixedParams.end(); ++itFix){
      upar.Fix( (*itFix) );
    }
    
    // 		bool prescan=false;
    // 		if(prescan){
    // 			upar.Fix(0);
    // 			MnScan theScan(theFcn, upar);
    // 			FunctionMinimum smin = theScan();
    // 			MnUserParameterState sState = smin.UserState();
    // 			cout << "After scan" << endl;
    // 			cout << sState << endl;
    
    // 			upar = smin.UserParameters();
    // 			upar.Release(0);
    // 		}
    
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
    fitParamsNew finalFitParams=theStartparams;
    theFitParamBase->getFitParamVal(finalParamVec, finalFitParams);
    
    // 		//MnUserCovariance theCov = min.UserCovariance() ;
    // 		//cout << "User vov : "<< endl;
    // 		//cout << theCov << endl;
    
    theFitParamBase->printParams(finalFitParams);
    double theLh=theLhPtr->calcLogLh(finalFitParams);
    Info <<"theLh = "<< theLh << endmsg;    
    
    const std::vector<double> finalParamErrorVec=finalUsrParameters.Errors();
    fitParamsNew finalFitErrs=theErrorparams;
    theFitParamBase->getFitParamVal(finalParamErrorVec, finalFitErrs);
    
    std::string finalResultname = "finalResult" + jobOption + ".dat";
    std::ofstream theStream ( finalResultname.c_str() );
    theFitParamBase->dumpParams(theStream, finalFitParams, finalFitErrs);

    
    MnUserCovariance theCovMatrix = min.UserCovariance();
    std::cout  << min << std::endl;
    
    //std::ofstream theCompStream ( "componentIntensity.dat");
    //theProdLh->dumpComponentIntensity( theCompStream, finalFitParams, theErrMatrix );
    JpsiGamEtaPiPiHistNew theHist(theLhPtr, finalFitParams,theAppParams.massRange(), jobOption.c_str());
    theHist.PrintToPDF(jobOption);
    theHist.SaveToROOT();
    end= clock();
    double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
    Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;




    // Start event number calculation for each wave
    std::cout << "Start event number calculation for each wave" << std::endl;    
   
    std::vector<std::string> hypVec_test=theAppParams.enabledHyps();
    fitParamsNew finalFitParams_test=finalFitParams;
    std::vector<std::string>::iterator it;
    int hypnumber=0;

    for (it=hypVec_test.begin(); it!=hypVec_test.end();++it){
      if (((*it).find("Gamma")==0 || (*it).find("Phasespace")==0)) hypnumber++;
    }
    std::cout << "Number of hypothesis found: " << hypnumber << std::endl;
    
    std::vector<double> evNumResult;
    std::string hypname;
    std::string suffix;

    int j;
    for (int i=1;i<=hypnumber;i++){
      j=1;
      for (it=hypVec_test.begin(); it!=hypVec_test.end();++it){
	// Mark bad hypothesis
        if (((*it).find("Gamma")==0 || (*it).find("Phasespace")==0 ) && i!=j) {
          (*it).insert(0, "#");
	  j++;
	}
	else{ 
	  if (((*it).find("Gamma")==0  || (*it).find("Phasespace")==0)){
	    j++;
	    //(*it).copy(hypname,(*it).size());
	    hypname = (*it).c_str();
	  }
	}
      }
      std::cout << "Start calulation with following hypothesis: " << hypname  << std::endl;
      if (startWithHyp=="production"){
	theLhPtr = boost::shared_ptr<AbsLhNew> (new JpsiGamEtaPiPiProdLhNew(theJpsiGamEtaPiPiEventListPtr, hypVec_test, jpsiGamEtaPiPiStatesPtr));
      }
      else {
	Alert << "start with hypothesis " << startWithHyp << " not supported!!!!" ;  // << endmsg;                                                                                                                                       
	exit(1);
      }
      suffix += jobOption.c_str();
      suffix += hypname.c_str();
      JpsiGamEtaPiPiHistNew theHist(theLhPtr, finalFitParams_test ,theAppParams.massRange(), suffix);
      suffix.clear();
      evNumResult.push_back(theHist.getFitEvents());
      theHist.SaveToROOT();
      // Unmark bad hypothesis
      for (it=hypVec_test.begin(); it!=hypVec_test.end();++it){
        if ((*it).find("#")==0) (*it).erase(0,1);
      }
      std::cout << std::endl;
    }
    


    // Global Summary Output
    int number_fitParams = upar.Params().size()-fixedParams.size();
    double a0_fitmass = 0;
    string aname;
    for (int i=0; i<int(upar.Params().size()); ++i){
      aname = upar.Name(i);
      if (!aname.find("a0_980Mass")) {
	std::cout << upar.Name(i) << "\t" << upar.Value(i) << "\t" << upar.Error(i) << std::endl;
	a0_fitmass =  upar.Value(i);
      }
    }

    std::ofstream summaryfile(sumFile.c_str(), std::ios::out|std::ios::app);
    summaryfile.precision(4);
    summaryfile << theAppParams.massRange().first  << "\t" << theAppParams.massRange().second << "\t" << theHist.getDataEvents() << "\t" << theHist.getMcEvents()  << "\t" << jobOption.c_str() << "\t" << theLh << "\t" << number_fitParams << "\t" << finalUsrParameters.Value("a0_980Mass");
    for (unsigned int i=0;i<evNumResult.size();i++){
      summaryfile << "\t" << evNumResult[i];
    }
    summaryfile <<  std::endl;
    summaryfile.close();


    return 0;
  }
  
  return 0;
}

