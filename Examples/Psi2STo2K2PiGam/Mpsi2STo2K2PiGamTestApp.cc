#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamParser.hh"

#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamHist.hh"
#include "Examples/Psi2STo2K2PiGam/AbsPsi2STo2K2PiGamLh.hh"
#include "Examples/Psi2STo2K2PiGam/HypProdLh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp1Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp2Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp3Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp4Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp5Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp6Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp7Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp8Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Hyp9Lh.hh"
#include "Examples/Psi2STo2K2PiGam/MPsi2STo2K2PiGamFcn.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "Examples/Psi2STo2K2PiGam/Stream2K2PiGamFitParms.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtGenConverter.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamHitAndMiss.hh"


#include "Setup/PwaEnv.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "Event/BesEvtReader.hh"
#include "Event/HepMCEventList.hh"
#include "Particle/PdtParser.hh"
#include "qft++/topincludes/tensor.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnPrint.h"

#include "PwaUtils/pbarpStates.hh"
#include "ErrLogger/ErrLogger.hh"

//#include "HepMC/GenEvent.h"


using namespace ROOT::Minuit2;

void setErrLogMode( const Psi2STo2K2PiGamParser::enErrLogMode& erlMode ) {
  switch(erlMode) {
  case Psi2STo2K2PiGamParser::debug :
    ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
    break;
  case Psi2STo2K2PiGamParser::trace :
    ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
    break;
  case Psi2STo2K2PiGamParser::routine :
    ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
    break;
  case Psi2STo2K2PiGamParser::warning :
    ErrLogger::instance()->setLevel(log4cpp::Priority::WARN);
    break;
  case Psi2STo2K2PiGamParser::error :
    ErrLogger::instance()->setLevel(log4cpp::Priority::ERROR);
    break;
  case Psi2STo2K2PiGamParser::alert :
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
  static Psi2STo2K2PiGamParser theAppParams(__argc, __argv);
  
  // Set the desired error logging mode
  setErrLogMode(theAppParams.getErrLogMode());
  
  // Set the Config File 
  std::string theCfgFile = theAppParams.getConfigFile();
  Info << "The path to config file is " << theCfgFile << "\n" << endmsg;    

  // Set the Fit Parameter file
  std::string paramStreamerPath=theAppParams.fitParamFile();
  Stream2K2PiGamFitParms theParamStreamer(paramStreamerPath);
  param2K2PiGam theStartparams=theParamStreamer.getFitParamVal();
  param2K2PiGam theErrorparams=theParamStreamer.getFitParamErr();

  bool genMode=theAppParams.generatorMode();

  EventList eventsData;
  EventList eventsMc;

  // Parsing events from data and mc file (if genMode = false)
  if (!genMode){
    const std::string datFile=theAppParams.dataFile();
    const std::string mcFile=theAppParams.mcFile();
    Info << "data file: " << datFile << endmsg;
    Info << "mc file: " << mcFile << endmsg;
    
    ParticleTable pTable;
    PdtParser parser;
    
    std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");
    std::string pdtFile(theSourcePath+"/Particle/pdt.table");
    if (!parser.parse(pdtFile, pTable)) {
      Alert << "Error: could not parse " << pdtFile << endmsg;
      exit(1);
    }
    
    std::vector<std::string> fileNames;
    
    fileNames.push_back(datFile);
    BesEvtReader eventReader(fileNames, 5, 0); 
    eventReader.fillAll(eventsData);
    
    if (!eventsData.findParticleTypes(pTable))
      Warning << "could not find all particles" << endmsg;
    Info << "\nFile has " << eventsData.size() << " events. Each event has "
	 <<  eventsData.nextEvent()->size() << " final state particles.\n" << endmsg;
    eventsData.rewind();
    
    std::vector<std::string> fileNamesMc;
    fileNamesMc.push_back(mcFile);
    BesEvtReader eventReaderMc(fileNamesMc, 5, 0); 
    eventReaderMc.fillAll(eventsMc);
    eventsMc.rewind();
  }

  boost::shared_ptr<const Psi2STo2K2PiGamEvtList> thePsi2STo2K2PiGamEvtListPtr(new Psi2STo2K2PiGamEvtList(eventsData, eventsMc));


  // Get enableHyp entries from Config File and disabling all hypotheses by default
  const std::vector<std::string> enabledHyps=theAppParams.enabledHyps();
  std::map<const std::string, bool> hypMap;
  hypMap["K1_1270Hyp"]=false;
  hypMap["K1_1400Hyp"]=false;
  hypMap["K0_1430_K0_1430Hyp"]=false;
  hypMap["K2_1430_K2_1430Hyp"]=false;
  hypMap["K0_1430_K2_1430Hyp"]=false;
  hypMap["K0_1430_K892Hyp1"]=false;
  hypMap["K2_1430_K892Hyp1"]=false;
  hypMap["K1_1410_K1_1410Hyp"]=false;
  hypMap["K1_1410_K892Hyp1"]=false;
  hypMap["f1710_f1710Hyp1"]=false;
  hypMap["doHyp2"]=false;
  hypMap["sigmaf980Hyp3"]=false;
  hypMap["sigmaf1710Hyp3"]=false;
  hypMap["sigmaf2200Hyp3"]=false;
  hypMap["f980f1370Hyp4"]=false;
  hypMap["f980f1500Hyp4"]=false;
  hypMap["f1710f1370Hyp4"]=false;
  hypMap["f980f_2_1270Hyp4"]=false;
  hypMap["f980f_2_1430Hyp4"]=false;
  hypMap["f980f_2_1525Hyp4"]=false;
  hypMap["f980f_2_1810Hyp4"]=false;
  hypMap["f980f_2_1950Hyp4"]=false;
  hypMap["f980f_2_2010Hyp4"]=false;
  hypMap["f1500f_2_1525Hyp4"]=false;
  hypMap["f1710f_2_1430Hyp4"]=false;
  hypMap["f1710f_2_1950Hyp4"]=false;
  hypMap["sigmaf1370Hyp4"]=false;
  hypMap["sigmaf_2_1270Hyp4"]=false;
  hypMap["sigmaf_2_1525Hyp4"]=false;
  hypMap["sigmaf_2_1810Hyp4"]=false;
  hypMap["sigmaf_2_1950Hyp4"]=false;
  hypMap["f_2_1270f_2_1270Hyp4"]=false;
  hypMap["f_2_1270f_2_1525Hyp4"]=false;
  hypMap["f_2_1270f_2_1810Hyp4"]=false;
  hypMap["f_2_1270f_2_1950Hyp4"]=false;
  hypMap["f_2_1810f_2_1810Hyp4"]=false;
  hypMap["f_2_1950f_2_1950Hyp4"]=false;
  hypMap["doHyp5"]=true;
  hypMap["K_0_2400KHyp5"]=false;
  hypMap["K_0_2400KTof_0_1710Hyp5"]=false;
  hypMap["K_1_2400KHyp5"]=false;
  hypMap["K_1_2400KTof_0_1710Hyp5"]=false;
  hypMap["K_1_2400KToK_0_1430Hyp5"]=false;
  hypMap["K_2_2400KTof980Hyp5"]=false;
  hypMap["K_2_2400KTof_0_1710Hyp5"]=false;
  hypMap["ChiToK_0_1430KPiHyp5"]=false;
  hypMap["ChiToK892KPiHyp5"]=false;
  hypMap["K_0_1430K_0_1950Hyp6"]=false;
  hypMap["KappaK892Hyp6"]=false;
  hypMap["KappaK_0_1430Hyp6"]=false;
  hypMap["KappaK_0_1950Hyp6"]=false;
  hypMap["K1_1680Hyp7"]=false;
  hypMap["K1_1680K1_1680Hyp7"]=false;
  hypMap["K1_1680K0_1430Hyp7"]=false;
  hypMap["K1_2300Hyp7"]=false; 
  hypMap["KappaHyp7"]=false;
  hypMap["K_0_1460ToKstPiHyp8"]=false;
  hypMap["K_0_1460ToSigmaKHyp8"]=false;
  hypMap["K_0_1830ToKstPiHyp8"]=false;
  hypMap["K_1_1650Hyp8"]=false;
  hypMap["Pi_2_1670Tof_2_1270PiHyp9"]=false;
  hypMap["Pi_2_1670ToK892KHyp9"]=false;
  hypMap["Pi1800Tof980PiHyp9"]=false;
  hypMap["Pi1800Tof1370PiHyp9"]=false;
  hypMap["Pi_2_2285Tof1710PiHyp9"]=false;
  hypMap["Pi_2_2285ToK892KHyp9"]=false;
  hypMap["Pi1800Pi0ToK892KHyp9"]=false;
  hypMap["Pi3000Pi0ToK892KHyp9"]=false;
  hypMap["Pi3000Pi0ToK_0_1950KHyp9"]=false;
  hypMap["Pi_2_2285ToK_0_1430KHyp9"]=false;
  hypMap["Pi_2_2285ToK_2_1430KHyp9"]=false;
  hypMap["f980f_2_2300Hyp9"]=false;
  hypMap["f_2_2300sigmaHyp9"]=false;
  hypMap["K_2_1770ToK_2_1430PiHyp9"]=false;
  hypMap["Pi1800ToKappaKHyp9"]=false;
  hypMap["K_2_1820ToK_2_1430PiHyp9"]=false;
  
  std::vector<std::string>::const_iterator itStr;
  

  // Enabling hypotheses due to enableHyp entries in Config File
  for (itStr=enabledHyps.begin(); itStr!=enabledHyps.end(); ++itStr){
    std::map<const std::string, bool>::const_iterator iter= hypMap.find( (*itStr) );
    if (iter !=hypMap.end()){
      hypMap[iter->first]= true;
      Info<< "hypothesis " << iter->first << " enabled" <<endmsg;
    }
    else { Alert << "hypothesis " << (*itStr) << " can not be enabled"<<endmsg;
      exit(0);
    }
  }

  boost::shared_ptr<Psi2STo2K2PiGamStates> theStatesPtr(new Psi2STo2K2PiGamStates());
  boost::shared_ptr<AbsPsi2STo2K2PiGamLh> thePsi2STo2K2PiGamLhPtr;


  // Choosing the hypothesis to start with
  std::string startWithHyp=theAppParams.startHypo();
  if (startWithHyp=="prod") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new HypProdLh(thePsi2STo2K2PiGamEvtListPtr, theStatesPtr));
  else if (startWithHyp=="hyp1") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp1Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap, theStatesPtr));
  else if (startWithHyp=="hyp2") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp2Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap, theStatesPtr));
  else if (startWithHyp=="hyp3") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp3Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap, theStatesPtr));
  else if (startWithHyp=="hyp4") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp4Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap, theStatesPtr));
  else if (startWithHyp=="hyp5") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp5Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap, theStatesPtr)); 
  else if (startWithHyp=="hyp6") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp6Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap, theStatesPtr));
  else if (startWithHyp=="hyp7") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp7Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap, theStatesPtr));
  else if (startWithHyp=="hyp8") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp8Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap, theStatesPtr));
  else if (startWithHyp=="hyp9") thePsi2STo2K2PiGamLhPtr= boost::shared_ptr<AbsPsi2STo2K2PiGamLh>(new Hyp9Lh(thePsi2STo2K2PiGamEvtListPtr, hypMap, theStatesPtr));
  else { Alert << "start with hypothesis " << startWithHyp << " not supported!!!!" << endmsg;
    exit(1);
  }


  // Generate Monte Carlo events on basis of phase-space distributed Monte Carlo events and a given parameter set
  if (genMode){
    const std::string hepMCinFile=theAppParams.hepMCinFile();
    const std::string addSuffix="f";
    const std::string hepMCoutFile=theAppParams.hepMCinFile()+addSuffix;
    Info << "HepMCin file: " << hepMCinFile << endmsg;
    Info << "HepMCout file: " << hepMCoutFile << endmsg;
    boost::shared_ptr<HepMCEventList> hepMCEvtListPtr(new HepMCEventList(hepMCinFile));

    hepMCEvtListPtr->rewind();

    Psi2STo2K2PiGamHitAndMiss theHitAndMiss(thePsi2STo2K2PiGamLhPtr, theStartparams, hepMCEvtListPtr, "Psi2STo2K2PiGamGen.root");     
 
    theHitAndMiss.dumpToHepMCAscii(hepMCoutFile);

    end= clock();
    double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
    Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;
    return 0;
  }


  // QA Mode to check the results with given start parameters _without_ minimizing
  bool qaMode=theAppParams.qaMode();
  std::cout << "qaMode: " << qaMode << std::endl;
  if (qaMode){
    thePsi2STo2K2PiGamLhPtr->printCurrentFitResult(theStartparams);
    double theLh=thePsi2STo2K2PiGamLhPtr->calcLogLh(theStartparams);
    Info << "theLh = " << theLh << endmsg;

    Psi2STo2K2PiGamHist Psi2STo2K2PiGamHist(thePsi2STo2K2PiGamLhPtr, theStartparams);

    end= clock();
    double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
    Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;
    return 0;
  }


  thePsi2STo2K2PiGamLhPtr->cacheAmplitudes(true);
  MnUserParameters upar;
  thePsi2STo2K2PiGamLhPtr->setMnUsrParams(upar, theStartparams, theErrorparams);


  // Fixing parameters due to mnParFix entries in Config File
  const std::vector<std::string> fixedParams=theAppParams.fixedParams();
  std::vector<std::string>::const_iterator itFix;
  for (itFix=fixedParams.begin(); itFix!=fixedParams.end(); ++itFix){
    upar.Fix( (*itFix) );
  }


  // Creating function to be minimized
  MPsi2STo2K2PiGamFcn mPsi2STo2K2PiGamFcn(thePsi2STo2K2PiGamLhPtr);


  // Scan Mode (varying one parameter with all others being fixed; x=paramValue, y=logLh)
  //   be sure to have a directory "scans/" for dumping the result files.  
  bool scanMode=theAppParams.scanMode();
  std::cout << "scanMode: " << scanMode << std::endl;
  if (scanMode){
    const std::vector<std::string> scanParams=theAppParams.scanParams();
    std::vector<std::string>::const_iterator itScan;
    std::vector<unsigned int> scanList;
    for (itScan=scanParams.begin(); itScan!=scanParams.end(); ++itScan){
      scanList.push_back(upar.Index(*itScan));
    }
  
    MnScan scan(mPsi2STo2K2PiGamFcn, upar, 1);
    cout << "Scan parameters: " << scan.Parameters() << endl;

    MnPlot plot;
    for(unsigned int i = 0; i < scanList.size(); i++) {
      cout << "Name of current parameter: " << upar.GetName(scanList[i]) << endl;
      char streamFile[400];
      strcpy(streamFile, "scans/");
      strcat(streamFile, upar.GetName(scanList[i]).c_str());
      ofstream scanstream(streamFile, std::ios::out | ios::trunc);
      std::vector<std::pair<double, double> > xy = scan.Scan(scanList[i], 100, 0., 0.);
      plot(xy);
      cout << "Scan finished for parameter: " << upar.GetName(scanList[i]) << endl;
      for(unsigned int i = 0; i<xy.size(); i++) {
	scanstream << xy[i].first << "\t" << xy[i].second << endl;
      } 
      scanstream.close();
    }
    std::cout<<scan.Parameters()<<std::endl;
    return 0;
  }

  
  // Calling Minimizer
  MnMigrad migrad(mPsi2STo2K2PiGamFcn, upar);
  Info <<"start migrad "<< endmsg;
  FunctionMinimum min = migrad();
  

  // Try another minimization method if no valid FunctionMinimum found 
  if(!min.IsValid()) {
      Info <<"FM is invalid, try with strategy = 2."<< endmsg;
      MnMigrad migrad2(mPsi2STo2K2PiGamFcn, min.UserState(), MnStrategy(2));
      min = migrad2();
  }
  
  // Retrieving fit results from minimizer
  MnUserParameters finalUsrParameters=min.UserParameters();
  const std::vector<double> finalParamVec=finalUsrParameters.Params();
  param2K2PiGam finalFitParams;
  thePsi2STo2K2PiGamLhPtr->setFitParamVal(finalFitParams, finalParamVec);

  thePsi2STo2K2PiGamLhPtr->cacheAmplitudes(false);


  // Creating histograms and NTuples with final results
  Psi2STo2K2PiGamHist Psi2STo2K2PiGamHist(thePsi2STo2K2PiGamLhPtr, finalFitParams);


  // Print final fit result and logLh
  thePsi2STo2K2PiGamLhPtr->printCurrentFitResult(finalFitParams);
  double theLh=thePsi2STo2K2PiGamLhPtr->calcLogLh(finalFitParams);
  Info <<"theLh = "<< theLh << endmsg;

  
  // Dumping final fit result to file
  const std::vector<double> finalParamErrorVec=finalUsrParameters.Errors();
  param2K2PiGam finalErrParams;
  thePsi2STo2K2PiGamLhPtr->setFitParamVal(finalErrParams, finalParamErrorVec);
  std::ofstream theStream ( "finalResult.dat");
  std::string theSuffix="Val"; 
  thePsi2STo2K2PiGamLhPtr->dumpCurrentResult(theStream, finalFitParams, theSuffix);
  theSuffix="Err"; 
  thePsi2STo2K2PiGamLhPtr->dumpCurrentResult(theStream, finalErrParams, theSuffix);
 

  // Printing CPU time
  end= clock();
  double cpuTime= (end-start)/ (CLOCKS_PER_SEC);
  Info << "cpuTime:\t" << cpuTime << "\tsec" << endmsg;
  return 0;
}

