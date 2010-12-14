#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEventList.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamHist.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamReader.hh"

#include "Examples/Psi2SToKpKmPiGam/AbsPsi2SToKpKmPiGamLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamBaseLh.hh"
#include "Examples/Psi2SToKpKmPiGam/PsiToChic1GamProdLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamHyp1Lh.hh"

#include "Examples/Psi2SToKpKmPiGam/MPsi2SToKpKmPiGamFcn.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"

// #include "Examples/Psi2SToKpKmPiGam/StreamFitParams.hh"
#include "Examples/Psi2SToKpKmPiGam/StreamKpKmPiGamFitParms.hh"

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

int main(int __argc,char *__argv[]){

  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
		    || strcmp( __argv[1], "--help" ) == 0 ) ){

    Info << "This application is a PWA fit for the reaction\n\n"
	 << "Psi(2S) -> Chi_c1 + gamma  ; Chi_c1 -> K+ K- pi0\n\n"
	 << "It makes use of the Helicity Formalism\n"
	 << "-K1400: spin hypothesis for K*1400 resonance (available so far: spin1, spin2, prod) (default spin1)\n" 
	 << "-msg <errorLogMode> you can choose the mode for the error logger\n\n"
	 << "-paramFile <path>\n"
	 << "i.e. with './Mpsi2SToKpKmPiGamApp -K1400 spin02 -paramFile ./fitParamSpin02.dat -msg debugging' \n"  
	 << endmsg;
    return 0;
  }

  
  std::string msgModeStr="default";
  std::string K1400SpinStr="spin1";
  std::string paramFilePathStr="default";
  std::string qaModeStr="no";
  std::string hypStr="base";

  // decode arguments
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-msg"){
      optind++;
      msgModeStr = __argv[optind];
      found=true;
    }
    if (sw=="-K1400"){
      optind++;
      K1400SpinStr = __argv[optind];
      found=true;
    }
    if (sw=="-paramFile"){
      optind++;
      paramFilePathStr = __argv[optind];
      found=true;
    }
    if (sw=="-qaMode"){
      optind++;
      qaModeStr = __argv[optind];
      found=true;
    }
    if (sw=="-hyp"){
      optind++;
      hypStr = __argv[optind];
      found=true;
    }
    if (!found){
      Warning << "Unknown switch: " 
            << __argv[optind] << endmsg;
      optind++;
    }
    
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
  }

  if(msgModeStr == "debugging") ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  else if(msgModeStr == "trace") ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
  else if(msgModeStr == "routine") ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
  else if(msgModeStr == "warning") ErrLogger::instance()->setLevel(log4cpp::Priority::WARN);
  else if(msgModeStr == "error")   ErrLogger::instance()->setLevel(log4cpp::Priority::ERROR); 
  else {
    ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
    Warning << "ErrorLogger not (properly) set -> Use mode 'DEBUG' " << endmsg;  
  }

  std::map<const std::string, bool> hypMap;
  hypMap["K0_1430HypBase"]=false;
  hypMap["K1_1410HypBase"]=false;
  hypMap["K2_1430HypBase"]=false;

  if (K1400SpinStr=="spin012"){
    hypMap["K0_1430HypBase"]=true;
    hypMap["K1_1410HypBase"]=true;
    hypMap["K2_1430HypBase"]=true;
  }
  else if (K1400SpinStr=="spin01"){
    hypMap["K0_1430HypBase"]=true;
    hypMap["K1_1410HypBase"]=true;
  }
  else if (K1400SpinStr=="spin02"){
    hypMap["K0_1430HypBase"]=true;
    hypMap["K2_1430HypBase"]=true;
  }
  else if (K1400SpinStr=="spin12"){
    hypMap["K1_1410HypBase"]=true;
    hypMap["K2_1430HypBase"]=true;
  }
  else if (K1400SpinStr=="spin0"){
    hypMap["K0_1430HypBase"]=true;
  }
  else if (K1400SpinStr=="spin1"){
    hypMap["K1_1410HypBase"]=true;
  }
  else if (K1400SpinStr=="spin2"){
    hypMap["K2_1430HypBase"]=true;
  }


  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");

  std::string paramStreamerPath=theSourcePath+"/Examples/Psi2SToKpKmPiGam/data/testStream.dat";

  if ( paramFilePathStr != "default") paramStreamerPath=paramFilePathStr;  

  StreamKpKmPiGamFitParms theParamStreamer(paramStreamerPath);
  paramKpKmPiGam theStartparams=theParamStreamer.getFitParamVal();
  paramKpKmPiGam theErrorparams=theParamStreamer.getFitParamErr();

  std::string datFile=theSourcePath+"/Examples/Psi2SToKpKmPiGam/data/fitvectorDATA_kkpi0.dat";
  std::string mcFile=theSourcePath+"/Examples/Psi2SToKpKmPiGam/data/fitvectorMC_chic1_kkpi0.dat"; 

  Info << "data file: " << datFile << endmsg;
  Info << "mc file: " << mcFile << endmsg;
  
  
  ParticleTable pTable;
  PdtParser parser;
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");
  if (!parser.parse(pdtFile, pTable)) {
    Alert << "Error: could not parse " << pdtFile << endmsg;
    exit(1);
  }
  

  std::vector<std::string> fileNames;

  fileNames.push_back(datFile);
  Psi2SToKpKmPiGamReader eventReader(fileNames, 4, 0); 
  EventList eventsData;
  eventReader.fillAll(eventsData);

  if (!eventsData.findParticleTypes(pTable))
    Warning << "could not find all particles" << endmsg;

  Info << "\nFile has " << eventsData.size() << " events. Each event has "
                  <<  eventsData.nextEvent()->size() << " final state particles.\n" << endmsg;
  eventsData.rewind();

  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = eventsData.nextEvent()) != 0 && evtCount < 20) {
    Info << "\n" 
                    << *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
                    << *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
                    << *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
                    << *(anEvent->p4(3)) << "\tm = " << anEvent->p4(3)->Mass() << "\n"
                    << endmsg;
    ++evtCount;
  }
  eventsData.rewind();


  std::vector<std::string> fileNamesMc;

  fileNamesMc.push_back(mcFile);
  Psi2SToKpKmPiGamReader eventReaderMc(fileNamesMc, 4, 0); 
  EventList eventsMc;
  eventReaderMc.fillAll(eventsMc);
  eventsMc.rewind();

  boost::shared_ptr<const Psi2SToKpKmPiGamEventList> thePsi2SToKpKmPiGamEventListPtr(new Psi2SToKpKmPiGamEventList(eventsData, eventsMc));


  boost::shared_ptr<AbsPsi2SToKpKmPiGamLh> thePsi2SToKpKmPiGamLhPtr;


  if (hypStr=="prod")  thePsi2SToKpKmPiGamLhPtr= boost::shared_ptr<AbsPsi2SToKpKmPiGamLh>(new PsiToChic1GamProdLh(thePsi2SToKpKmPiGamEventListPtr));
  else if (hypStr=="base") thePsi2SToKpKmPiGamLhPtr=boost::shared_ptr<AbsPsi2SToKpKmPiGamLh>(new Psi2SToKpKmPiGamBaseLh(thePsi2SToKpKmPiGamEventListPtr, hypMap));
  else if (hypStr=="hyp1") thePsi2SToKpKmPiGamLhPtr=boost::shared_ptr<AbsPsi2SToKpKmPiGamLh>(new Psi2SToKpKmPiGamHyp1Lh(thePsi2SToKpKmPiGamEventListPtr, hypMap));
  else {
    Alert << "hypothesis " << hypStr << " not supported!!!" << endmsg;
    exit(1);
  }

  std::cout << "qaModeStr: " << qaModeStr << std::endl;

  if (qaModeStr=="yes"){
    thePsi2SToKpKmPiGamLhPtr->printCurrentFitResult(theStartparams);
    double theLh=thePsi2SToKpKmPiGamLhPtr->calcLogLh(theStartparams);
    Info <<"theLh = "<< theLh << endmsg;

    Psi2SToKpKmPiGamHist thePsi2SToKpKmPiGamHist(thePsi2SToKpKmPiGamLhPtr, theStartparams);
    return 0;
  }

   MPsi2SToKpKmPiGamFcn mPsi2SToKpKmPiGamFcn(thePsi2SToKpKmPiGamLhPtr);

   MnUserParameters upar;
//   thePsi2SToKpKmPiGamLhPtr->setMnUsrParams(upar);
   thePsi2SToKpKmPiGamLhPtr->setMnUsrParams(upar, theStartparams, theErrorparams); 
   upar.Fix(1);
   upar.Fix(2);
// //   for (unsigned int i=0; i<5; i++) upar.Fix(i);

// //   unsigned int uparSize=upar.Params().size();
// //   for (unsigned int i=5; i< uparSize; i++) upar.Fix(i);;  

  MnMigrad migrad(mPsi2SToKpKmPiGamFcn, upar);
  
  Info <<"start migrad "<< endmsg;
  FunctionMinimum min = migrad();

  if(!min.IsValid()) {
    //try with higher strategy
    Info <<"FM is invalid, try with strategy = 2."<< endmsg;
    MnMigrad migrad2(mPsi2SToKpKmPiGamFcn, min.UserState(), MnStrategy(2));
    min = migrad2();
  }

  MnUserParameters finalUsrParameters=min.UserParameters();
  const std::vector<double> finalParamVec=finalUsrParameters.Params();
 
  paramKpKmPiGam finalFitParams;
  thePsi2SToKpKmPiGamLhPtr->setFitParamVal(finalFitParams, finalParamVec);


  
  Psi2SToKpKmPiGamHist thePsi2SToKpKmPiGamHist(thePsi2SToKpKmPiGamLhPtr, finalFitParams);

  thePsi2SToKpKmPiGamLhPtr->printCurrentFitResult(finalFitParams);
  double theLh=thePsi2SToKpKmPiGamLhPtr->calcLogLh(finalFitParams);
  Info <<"theLh = "<< theLh << endmsg;
  
  // print final fit result
  const std::vector<double> finalParamErrorVec=finalUsrParameters.Errors();
  for (size_t i=0; i<finalParamVec.size(); i++)
    {
      Info << "Value: " << finalParamVec[i] << "\t Error: " << finalParamErrorVec[i] << endmsg;
    }  


    return 0;
}

