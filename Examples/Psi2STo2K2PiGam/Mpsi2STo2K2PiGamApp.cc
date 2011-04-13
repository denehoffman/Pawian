#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamHist.hh"
//#include "Examples/Psi2STo2K2PiGam/GPODExpectationChecksT.hpp"

#include "Examples/Psi2STo2K2PiGam/AbsPsi2STo2K2PiGamLh.hh"
// #include "Examples/Psi2STo2K2PiGam/PsiToChic0GamProdLh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp1Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp2Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp3Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp4Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp5Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp6Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp7Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp8Lh.hh"

#include "Examples/Psi2STo2K2PiGam/MPsi2STo2K2PiGamFcn.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "Examples/Psi2STo2K2PiGam/Stream2K2PiGamFitParms.hh"
//#include "Examples/Psi2STo2K2PiGam/GPsi2STo2K2PiGamFit.hh"
//#include "Examples/Psi2STo2K2PiGam/GPsi2STo2K2PiGamIndividual.hh"

#include "Setup/PwaEnv.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "Event/BesEvtReader.hh"
// #include "Event/CBElsaReader.hh"
#include "Particle/PdtParser.hh"
#include "qft++/topincludes/tensor.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

#include "PwaUtils/pbarpStates.hh"

//#include "Minuit2/MnUserTransformation.h"
using namespace ROOT::Minuit2;

int main(int __argc,char *__argv[]){

  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
		    || strcmp( __argv[1], "--help" ) == 0 ) ){

    Info << "This application is a PWA fit for the reaction\n\n"
	 << "Psi(2S) -> Chi_c0 + gamma  ; Chi_c0 -> K+ K- pi0 pi0\n\n"
	 << "It makes use of the Helicity Formalism\n"
	 << "-K1400: spin hypothesis for K*1400 resonance (available so far: spin1, spin2, prod) (default spin1)\n" 
	 << "-msg <errorLogMode> you can choose the mode for the error logger\n\n"
	 << "-paramFile <path>\n"
	 << "i.e. with './Mpsi2STo2K2PiGamApp -K1400 spin02 -paramFile ./fitParamSpin02.dat -msg debugging' \n"  
	 << endmsg;
    return 0;
  }

  
  std::string msgModeStr="default";
  std::string paramFilePathStr="default";
  std::string qaModeStr="no";
  std::string K1400SpinStr="spin1";

  // decode arguments
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-msg"){
      optind++;
      msgModeStr = __argv[optind];
      found=true;
    }
    if (sw=="-paramFile"){
      optind++;
      paramFilePathStr = __argv[optind];
      found=true;
    }
    if (sw=="-K1400"){
      optind++;
      K1400SpinStr = __argv[optind];
      found=true;
    }
    if (sw=="-qaMode"){
      optind++;
      qaModeStr = __argv[optind];
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

  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");


   std::string paramStreamerPath=theSourcePath+"/Examples/Psi2STo2K2PiGam/data/startParamPsi2S.dat";
   if ( paramFilePathStr != "default") paramStreamerPath=paramFilePathStr;
  

   Stream2K2PiGamFitParms theParamStreamer(paramStreamerPath);
   param2K2PiGam theStartparams=theParamStreamer.getFitParamVal();
   param2K2PiGam theErrorparams=theParamStreamer.getFitParamErr();

//    std::string datFile=theSourcePath+"/Examples/Psi2STo2K2PiGam/data/chic0_data_pwasample.dat";
//    std::string mcFile=theSourcePath+"/Examples/Psi2STo2K2PiGam/data/chic0_signalmc_pwasample.dat"; 

//    std::string datFile=theSourcePath+"/Examples/Psi2STo2K2PiGam/data/101123_chic0_data_pwasample.dat";
//    std::string mcFile=theSourcePath+"/Examples/Psi2STo2K2PiGam/data/101123_chic0_signalmc_pwasample.dat";
   
   //std::string datFile=theSourcePath+"/Examples/Psi2STo2K2PiGam/data/101213_chic0_data_pwasample.dat";
   //std::string mcFile=theSourcePath+"/Examples/Psi2STo2K2PiGam/data/101213_chic0_signalmc_pwasample.dat";
   
   std::string datFile=theSourcePath+"/Examples/Psi2STo2K2PiGam/data/110310_chic0_data_pwasample.dat";
   std::string mcFile=theSourcePath+"/Examples/Psi2STo2K2PiGam/data/110310_chic0_signalmc_pwasample.dat";

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
  BesEvtReader eventReader(fileNames, 5, 0); 
  EventList eventsData;
  eventReader.fillAll(eventsData);

  if (!eventsData.findParticleTypes(pTable))
    Warning << "could not find all particles" << endmsg;

  Info << "\nFile has " << eventsData.size() << " events. Each event has "
                  <<  eventsData.nextEvent()->size() << " final state particles.\n" << endmsg;
  eventsData.rewind();


  std::vector<std::string> fileNamesMc;

  fileNamesMc.push_back(mcFile);
  BesEvtReader eventReaderMc(fileNamesMc, 5, 0); 
  EventList eventsMc;
  eventReaderMc.fillAll(eventsMc);
  eventsMc.rewind();

  boost::shared_ptr<const Psi2STo2K2PiGamEvtList> thePsi2STo2K2PiGamEvtListPtr(new Psi2STo2K2PiGamEvtList(eventsData, eventsMc));


  boost::shared_ptr<AbsPsi2STo2K2PiGamLh> thePsi2STo2K2PiGamLhPtr;


  std::map<const std::string, bool> hypMap;
  hypMap["K1_1270Hyp"]=true;
  hypMap["K1_1400Hyp"]=false;
  hypMap["K0_1430_K0_1430Hyp"]=true;
  hypMap["K2_1430_K2_1430Hyp"]=true;
  hypMap["K0_1430_K2_1430Hyp"]=true;
  hypMap["K1_1400_K1_1400Hyp"]=true;
  hypMap["disableHyp2"]=false;
  hypMap["sigmaf980Hyp3"]=true;
  hypMap["sigmaf1710Hyp3"]=true;
  hypMap["sigmaf2200Hyp3"]=true;
  hypMap["f980f1370Hyp4"]=true;
  hypMap["f980f1500Hyp4"]=true;
  hypMap["f1710f1370Hyp4"]=true;
  hypMap["f980f_2_1430Hyp4"]=true;
  hypMap["f980f_2_1525Hyp4"]=true;
  hypMap["f1500f_2_1525Hyp4"]=true;
  hypMap["f1710f_2_1430Hyp4"]=true;
  hypMap["disableHyp5"]=false;
  hypMap["K_0_1430K_0_1950Hyp6"]=true;
  hypMap["KappaK_0_1950Hyp6"]=false;
  hypMap["K1_1680Hyp7"]=true;
  hypMap["K1_2300Hyp7"]=true; 
  hypMap["KappaHyp7"]=true;
  hypMap["K_0_1460ToKstPiHyp8"]=true;
  hypMap["K_0_1830ToKstPiHyp8"]=true;
  hypMap["K_1_1650Hyp8"]=true;

//   if (K1400SpinStr=="prod")  thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new PsiToChic0GamProdLh(thePsi2STo2K2PiGamEvtListPtr));
  if (K1400SpinStr=="hyp1") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp1Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap));
  else if (K1400SpinStr=="hyp2") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp2Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap));
  else if (K1400SpinStr=="hyp3") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp3Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap));
  else if (K1400SpinStr=="hyp4") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp4Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap));
  else if (K1400SpinStr=="hyp5") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp5Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap)); 
  else if (K1400SpinStr=="hyp6") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp6Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap));
  else if (K1400SpinStr=="hyp7") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp7Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap));
  else if (K1400SpinStr=="hyp8") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp8Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap));
  else { Alert << "K1400 resonance with spin " << K1400SpinStr << " not supported!!!!" << endmsg;
    exit(1);
  }

  std::cout << "qaModeStr: " << qaModeStr << std::endl;

  if (qaModeStr=="yes"){
    thePsi2STo2K2PiGamLhPtr->printCurrentFitResult(theStartparams);
    double theLh=thePsi2STo2K2PiGamLhPtr->calcLogLh(theStartparams);
    Info <<"theLh = "<< theLh << endmsg;

    Psi2STo2K2PiGamHist Psi2STo2K2PiGamHist(thePsi2STo2K2PiGamLhPtr, theStartparams);
    return 0;
  }



  bool prefit=false;
  if (prefit){
//     MnUserParameters uparPre;
//     thePsi2STo2K2PiGamLhPtr->setMnUsrParams(uparPre, theStartparams, theErrorparams);
//     uparPre.Fix(1);
//     uparPre.Fix(2); 
//     GPsi2STo2K2PiGamFit preFit(thePsi2STo2K2PiGamLhPtr, uparPre);
//     preFit.doFit();
//     boost::shared_ptr<GPsi2STo2K2PiGamIndividual> bestIndividual=preFit.bestIndividual();
//     bestIndividual->getFitParams(theStartparams);
    Alert << "Geneva prefit is currently not implemented " << pdtFile << endmsg;
    exit(0);
  }

  MnUserParameters upar;
  thePsi2STo2K2PiGamLhPtr->setMnUsrParams(upar, theStartparams, theErrorparams);
  upar.Fix("J1P-1C-1L0S1ChiGamMag");
  upar.Fix("J1P-1C-1L0S1ChiGamphi");
  upar.Fix("J0P1C1L2S2K892K892Mag");
  upar.Fix("J0P1C1L2S2K892K892phi");
  
  if (hypMap["K1_1270Hyp"]==true){
    upar.Fix("J1P1C0L2S1K_1_1270ToK892PiMag");
    upar.Fix("J1P1C0L2S1K_1_1270ToK892Piphi");
  }
  
  if (hypMap["K2_1430_K2_1430Hyp"]==true){
    upar.Fix("J0P1C1L2S2K_2_1430K_2_1430Mag");
    upar.Fix("J0P1C1L2S2K_2_1430K_2_1430phi");
    upar.Fix("J0P1C1L4S4K_2_1430K_2_1430Mag");
    upar.Fix("J0P1C1L4S4K_2_1430K_2_1430phi");
  }

  if (hypMap["K1_1400_K1_1400Hyp"]==true){
    upar.Fix("J0P1C1L2S2K_1_1400K_1_1400Mag");
    upar.Fix("J0P1C1L2S2K_1_1400K_1_1400phi");
  }

  if (hypMap["K1_1680Hyp7"]==true){
    upar.Fix("J0P1C1L2S2K892K_1_1680Mag");
    upar.Fix("J0P1C1L2S2K892K_1_1680phi");
  }

  if (hypMap["K1_2300Hyp7"]==true){
    upar.Fix("J0P1C1L2S2K892K_1_2300Mag");
    upar.Fix("J0P1C1L2S2K892K_1_2300phi");
  }

  if (hypMap["K_1_1650Hyp8"]==true){
    upar.Fix("J1P1C0L2S1K_1_1650ToK892PiMag");
    upar.Fix("J1P1C0L2S1K_1_1650ToK892Piphi");
  }

  MPsi2STo2K2PiGamFcn mPsi2STo2K2PiGamFcn(thePsi2STo2K2PiGamLhPtr);
  
  MnMigrad migrad(mPsi2STo2K2PiGamFcn, upar);
  
  Info <<"start migrad "<< endmsg;
    FunctionMinimum min = migrad();
    
    if(!min.IsValid()) {
      //try with higher strategy
      Info <<"FM is invalid, try with strategy = 2."<< endmsg;
      MnMigrad migrad2(mPsi2STo2K2PiGamFcn, min.UserState(), MnStrategy(2));
      min = migrad2();
    }
    
    MnUserParameters finalUsrParameters=min.UserParameters();
    const std::vector<double> finalParamVec=finalUsrParameters.Params();
    
    param2K2PiGam finalFitParams;
    thePsi2STo2K2PiGamLhPtr->setFitParamVal(finalFitParams, finalParamVec);
    
    Psi2STo2K2PiGamHist Psi2STo2K2PiGamHist(thePsi2STo2K2PiGamLhPtr, finalFitParams);
    
    thePsi2STo2K2PiGamLhPtr->printCurrentFitResult(finalFitParams);
    double theLh=thePsi2STo2K2PiGamLhPtr->calcLogLh(finalFitParams);
    Info <<"theLh = "<< theLh << endmsg;
    
    // print final fit result
    const std::vector<double> finalParamErrorVec=finalUsrParameters.Errors();
    
    
    param2K2PiGam finalErrParams;
    thePsi2STo2K2PiGamLhPtr->setFitParamVal(finalErrParams, finalParamErrorVec);
    
    std::ofstream theStream ( "finalResult.dat");
    std::string theSuffix="Val"; 
    thePsi2STo2K2PiGamLhPtr->dumpCurrentResult(theStream, finalFitParams, theSuffix);
    theSuffix="Err"; 
    thePsi2STo2K2PiGamLhPtr->dumpCurrentResult(theStream, finalErrParams, theSuffix);
    
    return 0;
}

