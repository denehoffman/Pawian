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
#include "Examples/JpsiToPhiPhiGam/JpsiToPhiPhiGamEventList.hh"
#include "PwaUtils/PsiToXGamParser.hh"
#include "PwaUtils/PsiToXGamReader.hh"
#include "Examples/JpsiToPhiPhiGam/JpsiToPhiPhiGamHist.hh"
#include "Examples/JpsiToPhiPhiGam/JpsiToPhiPhiGamStates.hh"
//#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiEventListNew.hh"
#include "Examples/JpsiToPhiPhiGam/JpsiToPhiPhiGamLh.hh"

#include "PwaUtils/EvtDataBaseListNew.hh"
#include "PwaUtils/FitParamsBaseNew.hh"
#include "PwaUtils/StreamFitParmsBaseNew.hh"
#include "PwaUtils/PwaFcnBaseNew.hh"
#include "PwaUtils/AbsLhNew.hh"
#include "Utils/ErrLogUtils.hh"

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



int main(int __argc,char *__argv[]){
  clock_t start, end;
  start= clock();
  
  // Parse the command line
  static PsiToXGamParser theAppParams(__argc, __argv);
  
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
  
  std::pair<double, double> massRange = theAppParams.massRange();
  Info  << "Mass range: " << massRange.first << " " << massRange.second ;
  
  std::vector<std::string> fileNames;
  fileNames.push_back(datFile);

  bool withEvtWeight=theAppParams.useEvtWeight();
  Info << "EvtWeight: " << withEvtWeight << endmsg;
  PsiToXGamReader eventReader(fileNames, 5, 0, withEvtWeight);
  EventList eventsData;
  eventReader.fillMassRange(eventsData, massRange );
  
  if (!eventsData.findParticleTypes(pTable))
    Warning << "could not find all particles" ;  // << endmsg;
  
  Info 	<< "\nFile has " << eventsData.size() << " events. Each event has "
	<<  eventsData.nextEvent()->size() << " final state particles.\n" ;  // << endmsg;
  eventsData.rewind();
  
  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = eventsData.nextEvent()) != 0 && evtCount < 10) {
    Info	<< "\n"
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
  PsiToXGamReader eventReaderMc(fileNamesMc, 5, 0, withEvtWeight);
  EventList eventsMc;
  eventReaderMc.fillMassRange(eventsMc, massRange);
  eventsMc.rewind();
  
  //
  //calculate helicity angles, fill map with D-functions
  //
  
  boost::shared_ptr<JpsiToPhiPhiGamEventList> theJpsiGamXEventListPtr(new JpsiToPhiPhiGamEventList());
  theJpsiGamXEventListPtr->ratioMcToData(theAppParams.ratioMcToData());
  theJpsiGamXEventListPtr->read(eventsData, eventsMc);
  
  std::string mode=theAppParams.mode();
  std::cout << "Mode: " << mode << std::endl;
  if (mode=="plotmode"){
    // JpsiToPhiPhiGamHist theHist(theJpsiGamXEventListPtr,theAppParams.massRange()); 
    JpsiToPhiPhiGamHist theHist(theJpsiGamXEventListPtr);    
    // theHist.setMassRange(theAppParams.massRange() );
    return 0;
  }
  
  //
  //retrieve  hypotheses
  //
  
  boost::shared_ptr<JpsiToPhiPhiGamStates> jpsiToPhiPhiGamStatesPtr(new JpsiToPhiPhiGamStates()); 
  jpsiToPhiPhiGamStatesPtr->print(std::cout);

  const std::vector<std::string> hypVec=theAppParams.enabledHyps();
  boost::shared_ptr<AbsLhNew> theLhPtr;
  
  std::string startWithHyp=theAppParams.startHypo();
  
  if (startWithHyp=="production"){
    theLhPtr = boost::shared_ptr<AbsLhNew> (new JpsiToPhiPhiGamLh(theJpsiGamXEventListPtr, hypVec, jpsiToPhiPhiGamStatesPtr));
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
    
    theFitParamBase->dumpParams(theStreamDefault, defaultVal, defaultErr);
    return 0;
  }
  
  std::string paramStreamerPath=theAppParams.fitParamFile();
  StreamFitParmsBaseNew theParamStreamer(paramStreamerPath, theLhPtr);
  fitParamsNew theStartparams=theParamStreamer.getFitParamVal();
  fitParamsNew theErrorparams=theParamStreamer.getFitParamErr();

  PwaFcnBaseNew theFcn(theLhPtr, theFitParamBase);
  MnUserParameters upar;
  theFitParamBase->setMnUsrParams(upar, theStartparams, theErrorparams);
  
  std::cout << "\n\n**************** Minuit Fit parameter **************************" << std::endl;
  for (int i=0; i<int(upar.Params().size()); ++i){
    std::cout << upar.Name(i) << "\t" << upar.Value(i) << "\t" << upar.Error(i) << std::endl;
  } 

    const std::vector<std::string> fixedParams=theAppParams.fixedParams();  

    const unsigned int noOfFreeFitParams = upar.Params().size()-fixedParams.size();  
    if (mode=="qaMode"){

  Info << "\nThe parameter values are: " << "\n" << endmsg;
  theFitParamBase->printParams(theStartparams);
    
  Info << "\nThe parameter errors are: " << "\n" << endmsg;
  theFitParamBase->printParams(theErrorparams);
    
  double theLh=theLhPtr->calcLogLh(theStartparams);
  Info <<"theLh = "<< theLh << endmsg;
    
    
  JpsiToPhiPhiGamHist theHist(theLhPtr, theStartparams);
  theHist.setMassRange(theAppParams.massRange() );

  double evtWeightSumData = theJpsiGamXEventListPtr->NoOfWeightedDataEvts();
  double BICcriterion=2.*theLh+noOfFreeFitParams*log(evtWeightSumData);
  double AICcriterion=2.*theLh+2.*noOfFreeFitParams;
  double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)/(evtWeightSumData-noOfFreeFitParams-1);

  Info << "noOfFreeFitParams:\t" <<noOfFreeFitParams;
  Info << "evtWeightSumData:\t" <<evtWeightSumData; 
  Info << "BIC:\t" << BICcriterion << endmsg;
  Info << "AIC:\t" << AICcriterion << endmsg;
  Info << "AICc:\t" << AICccriterion << endmsg;

  std::string qaSummaryFileName = "qaSummary" + jobOption + ".dat";
  std::ofstream theQaStream ( qaSummaryFileName.c_str() );
  theQaStream << "BIC\t" << BICcriterion << "\n";
  theQaStream << "AICa\t" << AICcriterion << "\n";
  theQaStream << "AICc\t" << AICccriterion << "\n";
  theQaStream << "logLh\t" << theLh << "\n";
  theQaStream << "free parameter\t" << noOfFreeFitParams << "\n";
  theQaStream.close();
   
  end= clock();
  double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
  Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;
  return 0;
  }
  
  if (mode=="pwa"){
    bool cacheAmps = theAppParams.cacheAmps();
    Info << "caching amplitudes enabled / disabled:\t" <<  cacheAmps << endmsg;
    if (cacheAmps) theLhPtr->cacheAmplitudes();
    std::vector<std::string>::const_iterator itFix;
    for (itFix=fixedParams.begin(); itFix!=fixedParams.end(); ++itFix){
      upar.Fix( (*itFix) );
    }
    
  //   // 		bool prescan=false;
  //   // 		if(prescan){
  //   // 			upar.Fix(0);
  //   // 			MnScan theScan(theFcn, upar);
  //   // 			FunctionMinimum smin = theScan();
  //   // 			MnUserParameterState sState = smin.UserState();
  //   // 			cout << "After scan" << endl;
  //   // 			cout << sState << endl;
    
  //   // 			upar = smin.UserParameters();
  //   // 			upar.Release(0);
  //   // 		}
    
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
    //std::ofstream theStream ( "finalResult.dat");
    theFitParamBase->dumpParams(theStream, finalFitParams, finalFitErrs);
    
    MnUserCovariance theCovMatrix = min.UserCovariance();
    std::cout  << min << std::endl;
    
    // // 		//std::ofstream theCompStream ( "componentIntensity.dat");
    // // 		//theProdLh->dumpComponentIntensity( theCompStream, finalFitParams, theErrMatrix );
    // JpsiGamEtaPiPiHistNew theHist(theLhPtr, finalFitParams,theAppParams.massRange());
    // //theHist.setMassRange(theAppParams.massRange() );
    // // // 		theHist.fill();
    // end= clock();
    // double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
    // Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;




  //   std::cout << "Start event number calculation for each wave" << std::endl;    

  //   std::vector<std::string> hypVec_test=theAppParams.enabledHyps();
  //   fitParamsNew finalFitParams_test=finalFitParams;
  //   std::vector<std::string>::iterator it;
  //   int hypnumber=0;

  //   for (it=hypVec_test.begin(); it!=hypVec_test.end();++it){
  //     //std::cout << "hypothesis\t" << (*it) << "\t enabled\n";
  //     if ((*it).find("Gamma")==0) hypnumber++;
  //   }
  //   std::cout << std::endl;
  //   std::cout << "Number of hypothesis found: " << hypnumber << std::endl;

  //   int j;
  //   for (int i=1;i<=hypnumber;i++){
  //     j=1;
  //     for (it=hypVec_test.begin(); it!=hypVec_test.end();++it){
  // 	// Mark bad hypothesis
  //       if ((*it).find("Gamma")==0 && i!=j) {
  //         (*It).Insert(0, "#");
  // 	  j++;
  // 	}
  // 	else{ if ((*it).find("Gamma")==0) j++; }
  //     }
  //     std::cout << "Start calulation with following hypothesis:" << std::endl;
  //     if (startWithHyp=="production"){
  // 	theLhPtr = boost::shared_ptr<AbsLhNew> (new JpsiGamEtaPiPiProdLhNew(theJpsiGamEtaPiPiEventListPtr, hypVec_test, jpsiGamEtaPiPiStatesPtr));
  //     }
  //     else {
  // 	Alert << "start with hypothesis " << startWithHyp << " not supported!!!!" ;  // << endmsg;                                                                                                                                       
  // 	exit(1);
  //     }
  //     double theLh=theLhPtr->calcLogLh(finalFitParams_test);
  //     Info <<"theLh = "<< theLh << endmsg;
  //     JpsiGamEtaPiPiHistNew theHist(theLhPtr, finalFitParams_test ,theAppParams.massRange());
  //     // Unmark bad hypothesis
  //     for (it=hypVec_test.begin(); it!=hypVec_test.end();++it){
  //       if ((*it).find("#")==0) (*it).erase(0,1);
  //     }
  //     std::cout << std::endl;
  //   }
    
    return 0;
  }
  
  
  return 0;
}

