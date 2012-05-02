#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiParser.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiData.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiEventList.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiReader.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiHist.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiProdLh.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiFitParams.hh"

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
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnScan.h"

// #include "Examples/JpsiGamEtaPiPi/FitParamErrorMatrix.hh"
// #include "Examples/JpsiGamEtaPiPi/FitParamIndex.hh"
// #include "Examples/JpsiGamEtaPiPi/FitParamErrorMatrixStreamer.hh"



using namespace ROOT::Minuit2;


void setErrLogMode( const JpsiGamEtaPiPiParser::enErrLogMode& erlMode ) {
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
	clock_t start, end;
	start= clock();

	// Parse the command line
	static JpsiGamEtaPiPiParser theAppParams(__argc, __argv);

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
	boost::shared_ptr<const JpsiGamEtaPiPiEventList> theJpsiGamEtaPiPiEventListPtr(new JpsiGamEtaPiPiEventList(eventsData, eventsMc));
  
	//JpsiGamEtaPiPiHist theHist(theJpsiGamEtaPiPiEventListPtr);

	std::string mode=theAppParams.mode();
	std::cout << "Mode: " << mode << std::endl;
	if (mode=="plotmode"){
	  JpsiGamEtaPiPiHist theHist(theJpsiGamEtaPiPiEventListPtr);                                                               
	  theHist.setMassRange(theAppParams.massRange() );
	  return 0;
	}

	//
	//disable hypotheses
	//

	std::map<const std::string, bool> hypMap;
	hypMap["etaToPiPiEtaHyp"] =false;
	hypMap["etaToa980PiHyp"] =false;
	hypMap["etaToSigmaEtaHyp"] =false;
	hypMap["etaTof0_980EtaHyp"] =false;
	hypMap["etaToa2_1320PiHyp"] =false;
	hypMap["eta2ToPiPiEtaHyp"] =false;
	hypMap["eta2Toa980PiHyp"] =false;
	hypMap["eta2Toa2_1320PiHyp"] =false;
	hypMap["eta2ToSigmaEtaHyp"] =false;
	hypMap["eta2Tof0_980EtaHyp"] =false;
	hypMap["f1ToPiPiEtaHyp"] =false;
	hypMap["f1Toa980PiHyp"] =false;
	hypMap["usePhasespace"]=false;

	const std::vector<std::string> enabledHyps=theAppParams.enabledHyps();
	std::vector<std::string>::const_iterator itStr;

	for (itStr=enabledHyps.begin(); itStr!=enabledHyps.end(); ++itStr){
		std::map<const std::string, bool>::const_iterator iter= hypMap.find( (*itStr) );
		if (iter !=hypMap.end()){
			hypMap[iter->first]= true;
			Info<< "hypothesis " << iter->first << " enabled" ;  // << endmsg;
		}
		else { Alert << "hypothesis " << (*itStr) << " can not be enabled";  // << endmsg;
      	  	  exit(0);
		}
	}


	boost::shared_ptr<AbsLh> theLhPtr;

 	std::string startWithHyp=theAppParams.startHypo();

	if (startWithHyp=="production"){
		theLhPtr = boost::shared_ptr<AbsLh> (new JpsiGamEtaPiPiProdLh(theJpsiGamEtaPiPiEventListPtr, hypMap));
	}
	else {
		Alert << "start with hypothesis " << startWithHyp << " not supported!!!!" ;  // << endmsg;
		exit(1);
	}


	if (mode=="dumpDefaultParams"){
		fitParams defaultVal;
		fitParams defaultErr;
		theLhPtr->getDefaultParams(defaultVal, defaultErr);
		std::ofstream theStreamDefault ( "defaultparams.dat");
		boost::shared_ptr<FitParamsBase> theFitParamBase=boost::shared_ptr<FitParamsBase>(new JpsiGamEtaPiPiFitParams(defaultVal, defaultErr));
		theFitParamBase->dumpParams(theStreamDefault, defaultVal, defaultErr);
		return 0;
	}

	std::string paramStreamerPath=theAppParams.fitParamFile();

	StreamFitParmsBase theParamStreamer(paramStreamerPath, boost::shared_ptr<FitParamsBase> (new JpsiGamEtaPiPiFitParams()));
	fitParams theStartparams=theParamStreamer.getFitParamVal();
	fitParams theErrorparams=theParamStreamer.getFitParamErr();

// 	boost::shared_ptr<FitParamsBase> theFitParamBase
	boost::shared_ptr<FitParamsBase> theFitParamBase(new JpsiGamEtaPiPiFitParams(theStartparams, theErrorparams));

	if (mode=="qaMode"){

		Info << "\nThe parameter values are: " << "\n" << endmsg;
		theFitParamBase->printParams(theStartparams);

		Info << "\nThe parameter errors are: " << "\n" << endmsg;
		theFitParamBase->printParams(theErrorparams);

		double theLh=theLhPtr->calcLogLh(theStartparams);
		Info <<"theLh = "<< theLh << endmsg;

//        	std::string errFile = "finalErrorMatrix.dat";
// 		FitParamErrorMatrixStreamer theErrStreamer( errFile  );
// 		std::vector<double> theErrData;
// 		int ncols(0);
// 		theErrStreamer.matrixData( theErrData, ncols  );
// 		FitParamErrorMatrix theErrorMatrix(theErrData, ncols );
    
		JpsiGamEtaPiPiHist theHist(theLhPtr, theStartparams);
// 		theHist.setMassRange(theAppParams.massRange() );
// 		theHist.fill();

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

		bool prescan=false;
		if(prescan){
			upar.Fix(0);
			MnScan theScan(theFcn, upar);
			FunctionMinimum smin = theScan();
			MnUserParameterState sState = smin.UserState();
			cout << "After scan" << endl;
			cout << sState << endl;

			upar = smin.UserParameters();
			upar.Release(0);
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

		//MnUserCovariance theCov = min.UserCovariance() ;
		//cout << "User vov : "<< endl;
		//cout << theCov << endl;

		theFitParamBase->printParams(finalFitParams);
		double theLh=theLhPtr->calcLogLh(finalFitParams);
		Info <<"theLh = "<< theLh << endmsg;


		const std::vector<double> finalParamErrorVec=finalUsrParameters.Errors();
		fitParams finalFitErrs=theFitParamBase->getFitParamVal(finalParamErrorVec);

		std::ofstream theStream ( "finalResult.dat");
		theFitParamBase->dumpParams(theStream, finalFitParams, finalFitErrs);

		MnUserCovariance theCovMatrix = min.UserCovariance();
		std::cout  << min << std::endl;

		//std::ofstream theCompStream ( "componentIntensity.dat");
		//theProdLh->dumpComponentIntensity( theCompStream, finalFitParams, theErrMatrix );
 		JpsiGamEtaPiPiHist theHist(theLhPtr, finalFitParams);
// 		theHist.setMassRange(theAppParams.massRange() );
// 		theHist.fill();

		return 0;
	}


	return 0;
}

