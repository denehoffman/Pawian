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

std::map<std::string,double> doQaMode(fitParamsNew theStartparams, fitParamsNew theErrorparams,boost::shared_ptr<AbsLhNew> theLhPtr, std::pair<double, double> massRange, std::string jobOption, int noOfFreeFitParams){

  std::map<std::string,double> result;
  double theLh=theLhPtr->calcLogLh(theStartparams);
  Info <<"theLh = "<< theLh << endmsg;
  
  JpsiGamEtaPiPiHistNew theHist(theLhPtr, theStartparams, massRange, jobOption.c_str());
  theHist.PrintToPDF(jobOption);
  theHist.SaveToROOT();
    
  double evtWeightSumData = theHist.getDataEvents();
  double BICcriterion=2.*theLh+noOfFreeFitParams*log(evtWeightSumData);
  double AICcriterion=2.*theLh+2.*noOfFreeFitParams;
  double AICccriterion=AICcriterion+2.*noOfFreeFitParams*(noOfFreeFitParams+1)/(evtWeightSumData-noOfFreeFitParams-1);

  Info << "noOfFreeFitParams:\t" << noOfFreeFitParams;
  Info << "evtWeightSumData:\t" << evtWeightSumData;
  Info << "evtWeightSumFit:\t" << theHist.getFitEvents();
  Info << "BIC:\t" << BICcriterion << endmsg;
  Info << "AIC:\t" << AICcriterion << endmsg;
  Info << "AICc:\t" << AICccriterion << endmsg;

  result["theLh"] = theLh;
  result["dataEvents"] = evtWeightSumData;
  result["mcEvents"] = theHist.getMcEvents();
  result["fitEvents"] = theHist.getFitEvents();
  result["BIC"] = BICcriterion;
  result["AIC"] = AICcriterion;
  result["AICc"] = AICccriterion;

  return result;
}

void calcEventForEachWave(std::vector<std::string> hypVec, fitParamsNew finalFitParams, fitParamsNew finalFitErrs, boost::shared_ptr<JpsiGamEtaPiPiEventListNew> theJpsiGamEtaPiPiEventListPtr, boost::shared_ptr<JpsiGamEtaPiPiStates> jpsiGamEtaPiPiStatesPtr , std::pair<double, double> massRange, std::string jobOption, int noOfFreeFitParams, std::string sumFilePath, std::map<std::string,double> qaResult, double theLh){


  std::cout << "Start event number calculation for each wave" << std::endl;
  
  //  std::vector<std::string> hypVec=theAppParams.enabledHyps();
  std::vector<std::string>::iterator it;
  int hypnumber=0;
  
  // Count number of hypothesis                                                                                                                                                                                                            
  for (it=hypVec.begin(); it!=hypVec.end();++it){
    if (((*it).find("Gamma")==0 || (*it).find("Phasespace")==0)) hypnumber++;
  }
  std::cout << "Number of hypothesis found: " << hypnumber << std::endl;
  
  std::map<std::string,double> qaResultMod;
  std::vector<std::map<std::string,double> > allQaResults;
  std::vector<double> evNumResult;
  std::string hypname;
  std::string allhypname;
  std::string suffix;
  
  int j;
  // Manipulate hypothesis                                                                                                                                                                                                                 
  for (int i=1;i<=hypnumber;i++){
    j=1;
    for (it=hypVec.begin(); it!=hypVec.end();++it){
      // Mark bad hypothesis                                                                                                                                                                                                               
      if (((*it).find("Gamma")==0 || (*it).find("Phasespace")==0 ) && i!=j) {
	(*it).insert(0, "#");
	j++;
      }
      // Don't mark good hypothesis                                                                                                                                                                                                        
      else{
	if (((*it).find("Gamma")==0  || (*it).find("Phasespace")==0) && i==j){
	  j++;
	  //(*it).copy(hypname,(*it).size());                                                                                                                                                                                              
	  hypname = (*it).c_str();
	  allhypname+=hypname;
	}
      }
    }

    // Calculate with modified hypothesis                                                                                                                                                                                                  
    
    std::cout << "Start calulation with following hypothesis: " << hypname  << std::endl;
    boost::shared_ptr<AbsLhNew> theLhPtrWave;
    
    theLhPtrWave = boost::shared_ptr<AbsLhNew> (new JpsiGamEtaPiPiProdLhNew(theJpsiGamEtaPiPiEventListPtr, hypVec, jpsiGamEtaPiPiStatesPtr));
    
    suffix += jobOption.c_str();
    suffix += hypname.c_str();
    
    qaResultMod = doQaMode(finalFitParams, finalFitErrs, theLhPtrWave, massRange, suffix, noOfFreeFitParams);
    allQaResults.push_back(qaResultMod);
    
    std::cout << qaResult["fitEvents"] << "\t" << qaResult["dataEvents"] << std::endl;
    
    suffix.clear();
    
    // Unmark bad hypothesis                                                                                                                                                                                                               
    for (it=hypVec.begin(); it!=hypVec.end();++it){
      if ((*it).find("#")==0) (*it).erase(0,1);
    }
    std::cout << std::endl;
  }
  
  std::string sumFilename;
  int massRangemean = (int) ((massRange.first+(massRange.second-massRange.first)/2.)*1000);
  std::stringstream ss;
  ss << massRangemean;
  std::cout << massRangemean << std::endl;
  sumFilename+="summary";
  sumFilename+=ss.str();
  sumFilename+=".txt";
  std::string sumFile; sumFile+=sumFilePath; sumFile+=sumFilename;
  
  std::vector<std::map<std::string,double> >::iterator itsum;
  
  std::ofstream summaryfile(sumFile.c_str(), std::ios::out|std::ios::app);
  summaryfile.precision(4);
  summaryfile << massRangemean << "\t" << qaResult["dataEvents"] << "\t" << qaResult["mcEvents"] << "\t"  << qaResult["theLh"] << "\t" << noOfFreeFitParams << "\t"  << qaResult["AIC"] << "\t"  << qaResult["AICc"]<< "\t"  << qaResult["BIC"] << "\t" << allhypname << "\t" << jobOption.c_str();
  for (itsum=allQaResults.begin();itsum!=allQaResults.end() ;++itsum){
    summaryfile << "\t" << (*itsum)["fitEvents"];
  }
  
  summaryfile <<  std::endl;
  summaryfile.close();
  
  
}

int main(int __argc,char *__argv[]){
  //  setvbuf(stdout, NULL, _IONBF, 0); // Outputbuffer
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
  const std::string sumFilePath = theAppParams.sumFile();
  Info << "data file: " << datFile ;  // << endmsg;
  Info << "mc file: " << mcFile ;  // << endmsg;
  Info << "Summary file path: " << sumFilePath ;  // << endmsg;
  
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
  
  boost::shared_ptr<JpsiGamEtaPiPiEventListNew> theJpsiGamEtaPiPiEventListPtr(new JpsiGamEtaPiPiEventListNew());
  theJpsiGamEtaPiPiEventListPtr->ratioMcToData(theAppParams.ratioMcToData());
  theJpsiGamEtaPiPiEventListPtr->read(eventsData, eventsMc);
  
  std::string mode=theAppParams.mode();
  std::cout << "Mode: " << mode << std::endl;
  if (mode=="plotmode"){
    JpsiGamEtaPiPiHistNew theHist(theJpsiGamEtaPiPiEventListPtr,theAppParams.massRange(), jobOption.c_str());
    //    JpsiGamEtaPiPiHistNew theHist(theJpsiGamEtaPiPiEventListPtr);
    return 0;
  }
  
  //
  //retrieve  hypotheses
  //
  
  boost::shared_ptr<JpsiGamEtaPiPiStates> jpsiGamEtaPiPiStatesPtr(new JpsiGamEtaPiPiStates());
  jpsiGamEtaPiPiStatesPtr->print(std::cout);

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
    
    theFitParamBase->dumpParams(theStreamDefault, defaultVal, defaultErr);
    return 0;
  }
  
  std::string paramStreamerPath=theAppParams.fitParamFile();
  StreamFitParmsBaseNew theParamStreamer(paramStreamerPath, theLhPtr);
  fitParamsNew theStartparams=theParamStreamer.getFitParamVal();
  fitParamsNew theErrorparams=theParamStreamer.getFitParamErr();

  PwaFcnBaseNew theFcn(theLhPtr, theFitParamBase, jobOption);
  MnUserParameters upar;
  theFitParamBase->setMnUsrParams(upar, theStartparams, theErrorparams);
  
  std::cout << "\n\n**************** Minuit Fit parameter **************************" << std::endl;
  for (int i=0; i<int(upar.Params().size()); ++i){
    std::cout << upar.Name(i) << "\t" << upar.Value(i) << "\t" << upar.Error(i) << std::endl;
  }
  
  const std::vector<std::string> fixedParams=theAppParams.fixedParams();
  const unsigned int noOfFreeFitParams = upar.Params().size()-fixedParams.size(); 

  std::map<std::string,double> qaResult;
  
  if (mode=="qaMode"){
    Info << "\nThe parameter values are: " << "\n" << endmsg;
    theFitParamBase->printParams(theStartparams);
    Info << "\nThe parameter errors are: " << "\n" << endmsg;
    theFitParamBase->printParams(theErrorparams);

    qaResult = doQaMode(theStartparams, theErrorparams, theLhPtr, massRange, jobOption, noOfFreeFitParams);
    calcEventForEachWave(hypVec, theStartparams, theErrorparams, theJpsiGamEtaPiPiEventListPtr, jpsiGamEtaPiPiStatesPtr, massRange, jobOption, noOfFreeFitParams, sumFilePath, qaResult, qaResult["logLh"]);

    std::string qaSummaryFileName = "qaSummary" + jobOption + ".dat";
    std::ofstream theQaStream ( qaSummaryFileName.c_str() );
    theQaStream << "BIC\t" << qaResult["BIC"] << "\n";
    theQaStream << "AICa\t" << qaResult["AIC"] << "\n";
    theQaStream << "AICc\t" << qaResult["AICc"] << "\n";
    theQaStream << "logLh\t" << qaResult["logLh"] << "\n";
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

    qaResult = doQaMode(finalFitParams, finalFitErrs, theLhPtr, massRange, jobOption, noOfFreeFitParams);
    calcEventForEachWave(hypVec, finalFitParams, finalFitErrs, theJpsiGamEtaPiPiEventListPtr, jpsiGamEtaPiPiStatesPtr, massRange, jobOption, noOfFreeFitParams, sumFilePath, qaResult, theLh);
    
    end= clock();
    double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
    Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;

    
    ///////////////////////////////////////////////////
    // Start event number calculation for each wave  //
    ///////////////////////////////////////////////////

    // std::cout << "Start event number calculation for each wave" << std::endl;    
    
    // std::vector<std::string> hypVec_copy=theAppParams.enabledHyps();
    // std::vector<std::string>::iterator it;
    // int hypnumber=0;

    // // Count number of hypothesis
    // for (it=hypVec_copy.begin(); it!=hypVec_copy.end();++it){
    //   if (((*it).find("Gamma")==0 || (*it).find("Phasespace")==0)) hypnumber++;
    // }
    // std::cout << "Number of hypothesis found: " << hypnumber << std::endl;

    // std::map<std::string,double> qaResultMod;    
    // std::vector<std::map<std::string,double> > allQaResults;
    // std::vector<double> evNumResult;
    // std::string hypname;
    // std::string allhypname;
    // std::string suffix;
    
    // int j;
    // // Manipulate hypothesis
    // for (int i=1;i<=hypnumber;i++){
    //   j=1;
    //   for (it=hypVec_copy.begin(); it!=hypVec_copy.end();++it){
    //  	// Mark bad hypothesis
    // 	if (((*it).find("Gamma")==0 || (*it).find("Phasespace")==0 ) && i!=j) {
    // 	  (*it).insert(0, "#");
    //  	  j++;
    //  	}
    // 	// Don't mark good hypothesis
    //  	else{ 
    //  	  if (((*it).find("Gamma")==0  || (*it).find("Phasespace")==0) && i==j){
    //  	    j++;
    //  	    //(*it).copy(hypname,(*it).size());
    // 	    hypname = (*it).c_str();
    // 	    allhypname+=hypname;
    //  	  }
    //  	}
    //   }
      
    //   // Calculate with modified hypothesis

    //   std::cout << "Start calulation with following hypothesis: " << hypname  << std::endl;
    //   boost::shared_ptr<AbsLhNew> theLhPtrWave;
    //   if (startWithHyp=="production"){
    //  	theLhPtrWave = boost::shared_ptr<AbsLhNew> (new JpsiGamEtaPiPiProdLhNew(theJpsiGamEtaPiPiEventListPtr, hypVec_copy, jpsiGamEtaPiPiStatesPtr));
    //   }
    //   else {
    //  	Alert << "start with hypothesis " << startWithHyp << " not supported!!!!" ;  // << endmsg;                                                                                                                                       
    //  	exit(1);
    //   }
    //   suffix += jobOption.c_str();
    //   suffix += hypname.c_str();
      
    //   qaResultMod = doQaMode(finalFitParams, finalFitErrs, theLhPtrWave, massRange, suffix, noOfFreeFitParams);
    //   allQaResults.push_back(qaResultMod);

    //   std::cout << qaResult["fitEvents"] << "\t" << qaResult["dataEvents"] << std::endl;
            
    //   suffix.clear();
      
    //   // Unmark bad hypothesis
    //   for (it=hypVec_copy.begin(); it!=hypVec_copy.end();++it){
    // 	if ((*it).find("#")==0) (*it).erase(0,1);
    //   }
    //   std::cout << std::endl;
    // }
    
    // // Global Summary Output

    // //    int number_fitParams = upar.Params().size()-fixedParams.size();
    // //    double a0_fitmass = 0;
    // //    string aname;
    // //    for (int i=0; i<int(upar.Params().size()); ++i){
    // //      aname = upar.Name(i);
    // //      if (!aname.find("a0_980Mass")) {
    // //     	std::cout << upar.Name(i) << "\t" << upar.Value(i) << "\t" << upar.Error(i) << std::endl;
    // //     	a0_fitmass =  upar.Value(i);
    // //      }
    // //    }
    
    // std::string sumFilename;
    // int massRangemean = (int) ((massRange.first+(massRange.second-massRange.first)/2.)*1000);
    // std::stringstream ss;
    // ss << massRangemean;
    // std::cout << massRangemean << std::endl;
    // sumFilename+="summary";
    // sumFilename+=ss.str();
    // sumFilename+=".txt";
    // std::string sumFile; sumFile+=sumFilePath; sumFile+=sumFilename;

    // std::vector<std::map<std::string,double> >::iterator itsum;

    // std::ofstream summaryfile(sumFile.c_str(), std::ios::out|std::ios::app);
    // summaryfile.precision(4);
    // summaryfile << massRangemean << "\t" << allhypname << "\t" << jobOption.c_str() << "\t" << qaResult["dataEvents"] << "\t" << qaResult["mcEvents"] << "\t"  << theLh << "\t" << noOfFreeFitParams << "\t"  << qaResult["AIC"] << "\t"  << qaResult["AICc"]<< "\t"  << qaResult["BIC"];

    // for (itsum=allQaResults.begin();itsum!=allQaResults.end() ;++itsum){
    //   summaryfile << "\t" << (*itsum)["fitEvents"];
    // }
    // summaryfile <<  std::endl;
    // summaryfile.close();

    return 0;
  }
  
  return 0;
}

