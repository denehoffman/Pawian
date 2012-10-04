#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>

// Boost header files go here
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/logic/tribool.hpp"

#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"
#include "Examples/pbarpToOmegaPiLS/OmegaPiEventListLS.hh"
#include "Examples/pbarpToOmegaPiLS/OmegaTo3PiEventListLS.hh"
#include "Examples/pbarpToOmegaPiLS/OmegaPiHistLS.hh"
#include "Examples/pbarpToOmegaPiLS/pbarpToOmegaPi0StatesLS.hh"
#include "Examples/pbarpToOmegaPiLS/GArgumentParserLS.hh"
#include "Examples/pbarpToOmegaPiLS/MinuitstartparamLS.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiLhLS.hh"
#include "Examples/pbarpToOmegaPiLS/OmegaPiLhPi0GammaLS.hh"
#include "Examples/pbarpToOmegaPiLS/OmegaTo3PiLhPi0GammaLS.hh"
#include "Examples/pbarpToOmegaPiLS/OmegaTo3PiLhProdLS.hh"
#include "Examples/pbarpToOmegaPiLS/MOmegaPiFcnLS.hh"
#include "Examples/pbarpToOmegaPiLS/SpinDensityHistLS.hh"

// The individual that should be optimized
#include "Examples/pbarpToOmegaPiLS/GOmegaPiIndividualLS.hh"


#include "Setup/PwaEnv.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "Event/CBElsaReader.hh"
#include "Particle/PdtParser.hh"

#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/AbsStates.hh"
//#include "PwaUtils/pbarpStatesLS.hh"
#include "PwaUtils/DataUtils.hh"


#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

// GenEvA header files go here
#include <courtier/GAsioHelperFunctions.hpp>
#include <courtier/GAsioTCPClientT.hpp>
#include <courtier/GAsioTCPConsumerT.hpp>
#include <geneva/GBrokerEA.hpp>
#include <geneva/GEvolutionaryAlgorithm.hpp>
#include <geneva/GIndividual.hpp>
#include <geneva/GMultiThreadedEA.hpp>

/************************************************************************************************/
//Namespace aliases, so we do not need to quote the entire namespace name
namespace gg = Gem::Geneva;
namespace gp = Gem::Pawian;
namespace gc = Gem::Courtier;
namespace rm = ROOT::Minuit2;
namespace bp = boost::posix_time;
namespace bl = boost::logic;

/************************************************************************************************/
/**
 * This function constructs the path to the file.
 */
void constructPath(const std::string &thePrefix, const unsigned pbarMom, std::string &outFilePath)
{
  std::stringstream sstrDatFile; //String Stream for die construction of the path to the parameter File;
  sstrDatFile << thePrefix; 
  sstrDatFile.width(4);
  sstrDatFile.fill('0');
  sstrDatFile << right << pbarMom << ".dat";
  outFilePath = sstrDatFile.str();
}


/************************************************************************************************/
/**
 * This function checks if the file in the path theFilePath exists
 */
bool checkFileExist(const std::string &theFilePath)
{
  ifstream datChk(theFilePath.c_str());
  if (datChk) { return true; } 
  else { return false; }
}

//  
//  @param erlMode The desired error logging mode
//  

void setErrLogMode( const ApplicationParameterLS::enErrLogMode& erlMode ) {
  switch(erlMode) {
  case ApplicationParameterLS::debug :
    ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
    break;
  case ApplicationParameterLS::trace :
    ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
    break;
  case ApplicationParameterLS::routine :
    ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
    break;
  case ApplicationParameterLS::warning :
    ErrLogger::instance()->setLevel(log4cpp::Priority::WARN);
    break;
  case ApplicationParameterLS::error :
    ErrLogger::instance()->setLevel(log4cpp::Priority::ERROR);
    break;
  case ApplicationParameterLS::alert :
    ErrLogger::instance()->setLevel(log4cpp::Priority::ALERT);
    break;
  default: 
    ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  }
}

void emitExecutionMode( const ApplicationParameterLS::enExecMode& execMode ) {
  switch(execMode) {
    case ApplicationParameterLS::GenEvA:
      Info << "Using minimization algorithm: " << "GenEvA" << "\n" << endmsg;
      break;
    case ApplicationParameterLS::Minuit:
      Info << "Using minimization algorithm: " << "Minuit" << "\n" << endmsg;
      break;
    case ApplicationParameterLS::GenToMinuit:
      Info << "Using minimization algorithm: " << "First GenEvA then Minuit with final fit parameters from GenEvA" << "\n" << endmsg;
      break;
    case ApplicationParameterLS::SpinDensity:
      Info << "Calculating Spin Density" << "\n" << endmsg;
      break;
    case ApplicationParameterLS::QAmode:
      Info << "using QA mode" << "\n" << endmsg;
      break;
    case ApplicationParameterLS::StatesTest:
      Info << "using StatesTest mode" << "\n" << endmsg;
      break;
    case ApplicationParameterLS::HistTest:
      Info << "using HistTest mode" << "\n" << endmsg;
      break;
    case ApplicationParameterLS::ParamStreamTest:
      Info << "using test mode for streaming fit parameter" << "\n" << endmsg;
      break;
  }
}

bool streamFitParam(ApplicationParameterLS& appParam, MinuitStartParamLS& theUserPar){

  bool result=true;
  
  if ( !appParam.getPathStartParamFile().empty() )
    {
      ifstream theFile(appParam.getPathStartParamFile().c_str());   
      if (theFile)
        {
          Info << "Using start parameters from file " << appParam.getPathStartParamFile() << "\n" << endmsg;
          theUserPar.ParseStream(theFile);

	  Info << "The parameters are: " << endmsg;          
          theUserPar.printParamMap(std::cout);

        }
      else 
        {
          Warning << "Start parameter file " << appParam.getPathStartParamFile() << " not found!" << endmsg;
          result=false;
        }
    }
  else 
    {
      Info << "Start parameter file has to be set!" << endmsg;
      result=false;
    }

  return result;
}

bool streamFixMinuitParam(ApplicationParameterLS& appParam, std::vector<std::string>& paramNameVec){

  bool result=true;
  
  if ( !appParam.getMinuitFixParamFile().empty() )
    {
      ifstream theFile(appParam.getMinuitFixParamFile().c_str());   
      if (theFile)
        {
          Info << "Using file for fixing MINUIT paameters" << appParam.getMinuitFixParamFile() << "\n" << endmsg;

	  std::string strTmp;          
          
	  while (!theFile.eof()){
	    theFile >> strTmp;
	    if(theFile.eof()) break;
	    paramNameVec.push_back(strTmp); 
	  }

          Info << "The parameters fixed in MINUIT are as follows: " << "\n" << endmsg;

	  std::vector<std::string>::const_iterator it;
	  for (it=paramNameVec.begin(); it != paramNameVec.end(); ++it){
	    Info << (*it) << endmsg;
	  }
        }
      else 
        {
          Warning << "file for fixing MINUIT paameters " << appParam.getMinuitFixParamFile() << " not found!" << endmsg;
          result=false;
        }
    }
  else 
    {
      Info << "file for fixing MINUIT paameters has to be set!" << endmsg;
      result=false;
    }

  return result;
}

void states(boost::shared_ptr<const pbarpToOmegaPi0StatesLS> &theStates){
  theStates->print(std::cout);
}

void printEvents(EventList& theEvtList, unsigned int nEvtParticles, unsigned int nEvents){
  theEvtList.rewind();
  Event* anEvent;
  unsigned int evtCount = 0;
  while ((anEvent = theEvtList.nextEvent()) != 0 && evtCount < nEvents) {
    Info << "\n" << endmsg;
    for (unsigned int i=0; i<nEvtParticles;++i){
      Info << *(anEvent->p4(i)) << "\tm = " << anEvent->p4(i)->Mass() << "\n" << endmsg;
    }
    ++evtCount;
  }
  theEvtList.rewind();

}


void histTest(boost::shared_ptr<const AbsOmegaPiEventListLS>& theEvtList, const std::string &thePathToRootFile){
  OmegaPiHistLS(theEvtList, thePathToRootFile);
}

void paramStreamTest(ApplicationParameterLS& appParam, boost::shared_ptr<AbsOmegaPiLhLS> theOmLh){
  rm::MnUserParameters upar;
  MinuitStartParamLS theUserPar;

  bool setFitParam=streamFitParam(appParam, theUserPar);  
  if (!setFitParam) Warning << "using default start parameter" << endmsg;

  theOmLh->setMnUsrParams(upar, theUserPar);

  std::vector<double> theMinuitparams=upar.Params();
  
  Info << "\nminuit start parameters are:" << endmsg;

  for (unsigned int i=0; i<theMinuitparams.size(); ++i){
    std::cout << upar.GetName(i) << "\tvalue: " << upar.Value(i) << "\t error: " << upar.Error(i) << std::endl;   
  }

  OmegaPiDataLS::fitParamVal theParamVal;


  theOmLh->getFitParamVal(theParamVal, theMinuitparams);
  theOmLh->printFitParams(std::cout, theParamVal);

  Info << "dump the fit parameter if file streamTestFitStartParams.txt" << endmsg;
  std::ofstream theStream ("streamTestFitStartParams.txt");

  theOmLh->dumpCurrentResult(theStream, theParamVal);

  Info << "get again new MnUserParameters from current OmegaPiDataLS::fitParamVal" << endmsg;
  rm::MnUserParameters upar1; 

  theOmLh->setMnUsrParams(upar1, theParamVal);
  Info << "\nminuit parameters are:" << endmsg;

  std::cout << "upar1.Params().size()= " << upar1.Params().size() << std::endl;

  for (unsigned int i=0; i<upar1.Params().size(); ++i){
    std::cout << upar1.GetName(i) << "\tvalue: " << upar1.Value(i) << "\t error: " << upar1.Error(i) << std::endl;   
  }

  Info << "create default MnUserParameters" << endmsg;
  rm::MnUserParameters upar2; 
  theOmLh->setMnUsrParams(upar2);

  std::cout << "upar2.Params().size()= " << upar2.Params().size() << std::endl;

  for (unsigned int i=0; i<upar2.Params().size(); ++i){
    std::cout << upar2.GetName(i) << "\tvalue: " << upar2.Value(i) << "\t error: " << upar2.Error(i) << std::endl;   
  }  
}


bl::tribool Minuit(
     ApplicationParameterLS &theAppParams
     , OmegaPiDataLS::fitParamVal &fitParam
     , boost::shared_ptr<AbsOmegaPiLhLS> &omegaPiLh
     , bool streamParams
) {
  Info << "Minuit fit starts.\n" << endmsg;  

  // get pbarpToOmegaPi0States pointer back


  rm::MOmegaPiFcnLS mOmegaPiFcn(omegaPiLh);

  rm::MnUserParameters upar;

  if (streamParams){
    MinuitStartParamLS theUserPar;

    bool setFitParam=streamFitParam(theAppParams, theUserPar);  
    if (!setFitParam) Warning << "using default start parameter" << endmsg;
    
    omegaPiLh->setMnUsrParams(upar, theUserPar);
  }
  else{
    omegaPiLh->setMnUsrParams(upar, fitParam);
  }

  //now look if some parameters have to be fixed
  std::vector<std::string> fixedParams;

  if(streamFixMinuitParam(theAppParams, fixedParams)){
    omegaPiLh->fixMnUsrParams(upar, fixedParams);
  }


  rm::MnMigrad migrad(mOmegaPiFcn, upar);
  Info <<"start migrad "<< endmsg;
  rm::FunctionMinimum min = migrad();

  if(!min.IsValid()) {
    //try with higher strategy
    Info <<"FM is invalid, try with strategy = 2."<< endmsg;
    rm::MnMigrad migrad2(mOmegaPiFcn, min.UserState(), rm::MnStrategy(2));
    min = migrad2();
  }

  rm::MnUserParameters finalUsrParameters=min.UserParameters();

  // fill fitParam with the final fit result

  const std::vector<double> finalParamVec=finalUsrParameters.Params();
  omegaPiLh->getFitParamVal(fitParam, finalParamVec);
  cout << endl << "Final fit parameters:" << endl;
  omegaPiLh->printFitParams(cout, fitParam);

  cout << endl << "Errors:" << endl;
  const vector<double> finalParamErrors=finalUsrParameters.Errors();
  
  vector<double>::const_iterator it;
  int i;
  for (it = finalParamErrors.begin(), i=0; it < finalParamErrors.end(); it++,i++) 
    {
      std::string strName =finalUsrParameters.GetName(i);
      cout << "Name:" << strName << " ";
      cout << "Index:" << finalUsrParameters.Index(strName) << " ";
      cout << *it << endl;
    }

  return bl::tribool(true);
}







/*******************************************************************************************
*****/
/**
 * The actual optimization procedure using the Geneva library collection. Note that this fun
ction
 * will return a boost::logic::indeterminate value if this is a Geneva client in networked m
ode
 *
 * @return A boost::logic::tribool indicating whether the fit was successful or this is a cl
ient in networked mode
 */
bl::tribool GenEvA(
		   ApplicationParameterLS &theAppParams
		   , OmegaPiDataLS::fitParamVal &fitParm
		   , boost::shared_ptr<AbsOmegaPiLhLS> &omegaPiLh
		   ) {
  //--------------------------------------------------------------------------------------------
  // First check whether this is a client in networked mode. We then just start the listener and
  // return when it has finished.
  if(theAppParams.getParallelizationMode()==2 && !theAppParams.getServerMode()) {
    // Create a model that holds the static data needed for processing
    boost::shared_ptr<gp::GOmegaPiIndividualLS> gopi_ptr( new gp::GOmegaPiIndividualLS(omegaPiLh) );
//     boost::shared_ptr<gc::GAsioTCPClientT<gg::GIndividual> > p(
//        new gc::GAsioTCPClientT<gg::GIndividual>(
//              theAppParams.getIp()
//              , boost::lexical_cast<std::string>(theAppParams.getPort())
//           , gopi_ptr
//           )
//     );
//     p->setMaxStalls(0); // An infinite number of stalled data retrievals
//     p->setMaxConnectionAttempts(100); // Up to 100 failed connection attempts

//     // Prevent return of unsuccessful mutation attempts to the server
//     p->returnResultIfUnsuccessful(theAppParams.getReturnRegardless());

//     // Start the actual processing loop
//     p->run();

    return bl::indeterminate;
  }


  //--------------------------------------------------------------------------------------------
  Info << "GenEvA fit start.\n" << endmsg;

  // Create the first set of parent individuals. Initialization of parameters is done randomly.
  std::vector<boost::shared_ptr<gp::GOmegaPiIndividualLS> > parentIndividuals;
  for(std::size_t p = 0 ; p<theAppParams.getNParents(); p++) {
    boost::shared_ptr<gp::GOmegaPiIndividualLS> gopi_ptr( new gp::GOmegaPiIndividualLS(omegaPiLh) );
    gopi_ptr->setProcessingCycles(theAppParams.getProcessingCycles());
    // Make sure we start with a unique individual
    gopi_ptr->randomInit();
    
    parentIndividuals.push_back(gopi_ptr);
  }

  //----------------------------------------------------------------------------------------------
  // We can now start creating populations. We refer to them through the base class

  // This smart pointer will hold the different population types
  boost::shared_ptr<GEvolutionaryAlgorithm> pop_ptr;

  // Create the actual populations
  switch (theAppParams.getParallelizationMode()) {
  case 0: // Serial execution
    {
      // Create an empty population
      pop_ptr = boost::shared_ptr<GEvolutionaryAlgorithm>(new GEvolutionaryAlgorithm());
    }
    break;


  case 1: // Multi-threaded execution
    {
      // Create the multi-threaded population
      boost::shared_ptr<GMultiThreadedEA> popPar_ptr(new GMultiThreadedEA());

      // Population-specific settings
      popPar_ptr->setNThreads(theAppParams.getNEvaluationThreads());

      // Assignment to the base pointer
      pop_ptr = popPar_ptr;
    }
    break;

  case 2: // Networked execution (server-side)
    {
      // Create a network consumer and enrol it with the broker
      boost::shared_ptr<gc::GAsioTCPConsumerT<gg::GIndividual> > gatc(new gc::GAsioTCPConsumerT<gg::GIndividual>(theAppParams.getPort()));
      gatc->setSerializationMode(theAppParams.getSerMode());
      GINDIVIDUALBROKER->enrol(gatc);

      // Create the actual broker population
      boost::shared_ptr<GBrokerEA> popBroker_ptr(new GBrokerEA());
      popBroker_ptr->setWaitFactor(theAppParams.getWaitFactor());

      // Assignment to the base pointer
      pop_ptr = popBroker_ptr;
    }
    break;
  }

  //----------------------------------------------------------------------------------------------
  // Now we have suitable populations and can fill them with data

  // Add individuals to the population
  for(std::size_t p = 0 ; p<theAppParams.getNParents(); p++) {
    pop_ptr->push_back(parentIndividuals[p]);
  }

  // Specify some general population settings
  pop_ptr->setDefaultPopulationSize(theAppParams.getPopulationSize(),theAppParams.getNParents());
  pop_ptr->setMaxIteration(theAppParams.getMaxIterations());
  pop_ptr->setMaxTime(boost::posix_time::minutes(theAppParams.getMaxMinutes()));
  pop_ptr->setReportIteration(theAppParams.getReportIteration());
  pop_ptr->setRecombinationMethod(theAppParams.getRScheme());
  pop_ptr->setSortingScheme(theAppParams.getSmode());
  
  // Do the actual optimization
  pop_ptr->optimize();

  for(std::size_t i=0; i<pop_ptr->size(); i++) {
    std::cout << i << ": " << pop_ptr->at(i)->fitness() << std::endl;
  }

  //----------------------------------------------------------------------------------------------

  boost::shared_ptr<gp::GOmegaPiIndividualLS> bestIndividual_ptr=pop_ptr->getBestIndividual<gp::GOmegaPiIndividualLS>();

  bool fitParamsSet=bestIndividual_ptr->getFitParams(fitParm);
  
  if (!fitParamsSet){
    std::cout << "fit parameters could not set properly " << std::endl;
    exit(0);
  }
  
  
  omegaPiLh=bestIndividual_ptr->omegaPiLhPtr();

  Info << "GenEvA done.\n" << endmsg;

  return bl::tribool(true);

}

bl::tribool QAmode(ApplicationParameterLS &theAppParams,
                   OmegaPiDataLS::fitParamVal &finalFitParm,
                   boost::shared_ptr<AbsOmegaPiLhLS> &finalOmegaPiLh
                   )
{
  Info << "QA mode start.\n" << endmsg;
    
  rm::MOmegaPiFcnLS mOmegaPiFcn(finalOmegaPiLh);

  rm::MnUserParameters upar;

  MinuitStartParamLS theUserPar;

  bool setFitParam=streamFitParam(theAppParams, theUserPar);  
  if (!setFitParam){
    Alert << "parameter file not found" << endmsg;
    exit(0);
  }
  
  finalOmegaPiLh->setMnUsrParams(upar, theUserPar);
 
  const std::vector<double> finalParamVec=upar.Params();
  finalOmegaPiLh->getFitParamVal(finalFitParm, finalParamVec);
  cout << endl << " fit parameters:" << endl;
  finalOmegaPiLh->printFitParams(cout, finalFitParm);

  return bl::tribool(true);
}


bl::tribool calcSpinDensity(ApplicationParameterLS &theAppParams, 
                     boost::shared_ptr<AbsOmegaPiLhLS> finalOmegaPiLh,
                     OmegaPiDataLS::fitParamVal &finalFitParm)
{
  Info << "Starting spin density calculation." << "\n" << endmsg;
  OmegaPiDataLS::fitParamVal theParamVal;
        
  rm::MOmegaPiFcnLS mOmegaPiFcn(finalOmegaPiLh);
  rm::MnUserParameters upar;

  MinuitStartParamLS theUserPar;

  bool setFitParam=streamFitParam(theAppParams, theUserPar);  
  if (!setFitParam){
    Alert << "parameter file not found" << endmsg;
    exit(0);
  }

  finalOmegaPiLh->setMnUsrParams(upar, theUserPar);
  finalOmegaPiLh->getFitParamVal(theParamVal, upar.Params());       
  
  Info << "Using following fit parameter:\n" << endmsg;
  finalOmegaPiLh->printFitParams(std::cout, theParamVal);

  std::ostringstream theSpinDensityRootFile;

  if (theAppParams.getCalcAllSpindensity())
    {
      Info << "Calculating all spin density elements.\n" << endmsg;
      theSpinDensityRootFile << "./" << theAppParams.getName() << "SpinDensity" << "_Lmax" << theAppParams.getLMax() << "_mom" << theAppParams.getPbarMom() << ".root";

      SpinDensityHistLS theSpinDensityHist(theSpinDensityRootFile.str(), finalOmegaPiLh, theParamVal);
      theSpinDensityHist.createHistograms();
    }
  else
    {
      Info << "Calculating spin density elements for M=" << theAppParams.getM() << " M'=" << theAppParams.getM_() << ".\n" << endmsg;
      theSpinDensityRootFile << "./" << theAppParams.getName() << "SpinDensity_M1" << theAppParams.getM() << "_M2" << theAppParams.getM_() << "_Lmax" << theAppParams.getLMax() << "_mom" << theAppParams.getPbarMom() << ".root";
       
      SpinDensityHistLS theSpinDensityHist(theSpinDensityRootFile.str(), finalOmegaPiLh, theParamVal);
      theSpinDensityHist.createHistogram(theAppParams.getM(),theAppParams.getM_());
    }
  Info << "Spin density calculation done." << "\n" << endmsg;

  return bl::tribool(false);
}

int main(int __argc,char *__argv[]){

  // Disable output buffering
  setvbuf(stdout, NULL, _IONBF, 0);

  // Parse the command line
  static ApplicationParameterLS theAppParams(__argc, __argv);
  // Set the desired error logging mode
  setErrLogMode(theAppParams.getErrLogMode());

  // Inform the audience about the execution mode
  emitExecutionMode(theAppParams.getAppExecMode());

  std::string theCfgFile = theAppParams.getConfigFile();
  Info << "The path to config file is " << theCfgFile << "\n" << endmsg;
  Info << "Maximum orbital momentum content: " << theAppParams.getLMax() << endmsg;
  Info << "pbar momentum: " << theAppParams.getPbarMom()   << endmsg;

  std::string theFixParamMinuitFile = theAppParams.getMinuitFixParamFile();
  Info << "\n" << "The path to file for fixing MINUIT parameters is " << theFixParamMinuitFile << "\n" << endmsg;


  int lmax=theAppParams.getLMax();
  boost::shared_ptr<const pbarpToOmegaPi0StatesLS> pbarpToOmegaPi0StatesPtr(new pbarpToOmegaPi0StatesLS(lmax));

  if (theAppParams.getAppExecMode()==ApplicationParameterLS::StatesTest){
     states(pbarpToOmegaPi0StatesPtr);
     exit(1);
   }

  std::string piomegaDatFile=theAppParams.dataFile();
  std::string piomegaMcFile=theAppParams.mcFile(); 
  Info << "data file: " << piomegaDatFile << endmsg;
  Info << "mc file: " << piomegaMcFile << endmsg;
  
  int nParticlesPerEvt=0;
  bool readWeightData=true;

  if (theAppParams.getLhMode()=="OmegaPiLhGamma" || (theAppParams.getLhMode()=="OmegaPiLhGammaBw") ){
     nParticlesPerEvt=3;
   }
   else if (theAppParams.getLhMode()=="OmegaTo3PiLhGamma" || (theAppParams.getLhMode()=="OmegaTo3PiLhProd") || (theAppParams.getLhMode()=="OmegaTo3PiLhGammaBw")){
     nParticlesPerEvt=4;
   }

  else {
    Alert << "LhMode: " << theAppParams.getLhMode() << " does not exist!!!" << endmsg;
    exit(0);
  }

  
  if (checkFileExist(piomegaDatFile)){
      Info << "Using Data file " << piomegaDatFile << endmsg;
  }
  else{
      Alert <<"Data file for pbarMom= " << theAppParams.getPbarMom() << " not available !"  << endmsg;
      Alert << "File " << piomegaDatFile << " is missing !" << endmsg;
      exit(1);
  }

  if (checkFileExist(piomegaMcFile)){
      Info << "Using Monte Carlo file " << piomegaMcFile << endmsg;
  }
  else{
      Alert <<"Monte Carlo file for pbarMom= " << theAppParams.getPbarMom() << "not available !" << endmsg;
      Alert << "File " << piomegaMcFile << " is missing !" << endmsg;
      exit(1);
  }
  
  Info << "data file: " << piomegaDatFile << endmsg;
  Info << "mc file: " << piomegaMcFile << endmsg;

  ParticleTable pTable;
  PdtParser parser;

  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");
  if (!parser.parse(pdtFile, pTable)) {
    Alert << "Error: could not parse " << pdtFile << endmsg;
    exit(1);
  }

  std::vector<std::string> fileNames;
  fileNames.push_back(piomegaDatFile);
  EventList theDataEventList;

  boost::shared_ptr<CBElsaReader> eventReaderPtr;

  std::vector<std::string> fileNamesMc;

  fileNamesMc.push_back(piomegaMcFile);
  boost::shared_ptr<CBElsaReader> eventReaderMCPtr;
  EventList theMcEventList;

  eventReaderPtr = boost::shared_ptr<CBElsaReader>(new CBElsaReader(fileNames,nParticlesPerEvt, 0, readWeightData));
  eventReaderPtr->fillAll(theDataEventList);

  Info << "\nFile has " << theDataEventList.size() << " events. Each event has "
       <<  theDataEventList.nextEvent()->size() << " final state particles.\n" << endmsg;

  if (!theDataEventList.findParticleTypes(pTable)) {
    Warning << "could not find all particles" << endmsg;
  }
  theDataEventList.rewind();

  Info << "\nThe first data events are: " << endmsg;
  printEvents(theDataEventList,nParticlesPerEvt , 20);


  eventReaderMCPtr= boost::shared_ptr<CBElsaReader>(new CBElsaReader(fileNamesMc, nParticlesPerEvt, 0));
  eventReaderMCPtr->fillAll(theMcEventList);
  theMcEventList.rewind();

  Info << "\nThe first MC events are: " << endmsg;
  printEvents(theMcEventList,nParticlesPerEvt , 20);


  std::ostringstream theRootFilePath;
  boost::shared_ptr<const AbsOmegaPiEventListLS> theOmegaPiEventPtr;

  if (theAppParams.getLhMode()=="OmegaPiLhGamma" || (theAppParams.getLhMode()=="OmegaPiLhGammaBw") ){
    theOmegaPiEventPtr = boost::shared_ptr<const AbsOmegaPiEventListLS>(new OmegaPiEventListLS (theDataEventList, theMcEventList, theAppParams.getLMax()+1,  theAppParams.getPbarMom() ) ); 
    theRootFilePath << "./" << theAppParams.getName() << "OmegaPi0Fit_Lmax" << theAppParams.getLMax() << "_mom" << theAppParams.getPbarMom() << ".root";
  }
  else if ( (theAppParams.getLhMode()=="OmegaTo3PiLhGamma") || (theAppParams.getLhMode()=="OmegaTo3PiLhProd") || (theAppParams.getLhMode()=="OmegaTo3PiLhGammaBw") ){
    theOmegaPiEventPtr = boost::shared_ptr<const AbsOmegaPiEventListLS>(new OmegaTo3PiEventListLS (theDataEventList, theMcEventList, theAppParams.getLMax()+1,  theAppParams.getPbarMom() ) );
   theRootFilePath << "./" << theAppParams.getName() << "OmegaTo3PiFit_Lmax" << theAppParams.getLMax() << "_mom" << theAppParams.getPbarMom() << ".root";
  }
  if (theAppParams.getAppExecMode()==ApplicationParameterLS::HistTest){
    histTest(theOmegaPiEventPtr, theRootFilePath.str());
    exit(1);
  }

//   rm::MnUserParameters upar;
//   MinuitStartParamLS theUserPar;

//   bool setFitParam=streamFitParam(theAppParams, theUserPar); 
  


  boost::shared_ptr<AbsOmegaPiLhLS> theOmegaPiLh;
  if (theAppParams.getLhMode()=="OmegaPiLhGamma"){
    theOmegaPiLh = boost::shared_ptr<AbsOmegaPiLhLS>(new OmegaPiLhPi0GammaLS(theOmegaPiEventPtr, pbarpToOmegaPi0StatesPtr));  
  }
  else if (theAppParams.getLhMode()=="OmegaPiLhGammaBw"){
//     theOmegaPiLh = boost::shared_ptr<AbsOmegaPiLhLS>(new OmegaPiLhPi0GammaLS(theOmegaPiEventPtr, pbarpToOmegaPi0StatesPtr));
    Alert << "LhMode OmegaPiLhGammaBw not implemented so far!!!" << endmsg; 
  }
  else if (theAppParams.getLhMode()=="OmegaTo3PiLhGamma"){
    theOmegaPiLh = boost::shared_ptr<AbsOmegaPiLhLS>(new OmegaTo3PiLhPi0GammaLS(theOmegaPiEventPtr, pbarpToOmegaPi0StatesPtr));
  }
  else if (theAppParams.getLhMode()=="OmegaTo3PiLhProd"){
    theOmegaPiLh = boost::shared_ptr<AbsOmegaPiLhLS>(new OmegaTo3PiLhProdLS(theOmegaPiEventPtr, pbarpToOmegaPi0StatesPtr));
  }
  else{
    Alert <<"LhMode " << theAppParams.getLhMode() << " doesn't exist !!! " << endmsg;
    exit(1);
  }

  if (theAppParams.getAppExecMode()==ApplicationParameterLS::ParamStreamTest){
    paramStreamTest(theAppParams, theOmegaPiLh);
    exit(1);
  }

  bl::tribool bExecFinish=bl::tribool(false);
  OmegaPiDataLS::fitParamVal theParamVal;

  switch(theAppParams.getAppExecMode()) {  

  case ApplicationParameterLS::Minuit:
    {
      bExecFinish = Minuit(theAppParams, theParamVal, theOmegaPiLh, true);

      //now dump the final fit result
      std::ostringstream theParamFilePathMin;
      theParamFilePathMin << "./" << theAppParams.getName() << "LastFitParamOmegaPi0Fit_Mmax" << theAppParams.getLMax() << "_mom" << theAppParams.getPbarMom() << ".txt";
      ofstream outfileMin (theParamFilePathMin.str().c_str());
      theOmegaPiLh->dumpCurrentResult(outfileMin, theParamVal);
      outfileMin.close();
      break;
    }

  case ApplicationParameterLS::GenEvA:
    {
      bExecFinish = GenEvA(theAppParams, theParamVal, theOmegaPiLh);
      if(bl::indeterminate(bExecFinish)) return 0; // Simply terminate -- this is a client in networked mode that has done its job
    }

  case ApplicationParameterLS::GenToMinuit:
    {
      bExecFinish = GenEvA(theAppParams, theParamVal, theOmegaPiLh);
      Info << "Obtained NLL with GenEvA:\n" << theOmegaPiLh->calcLogLh(theParamVal) << endmsg;
      Info << "\n and fit parameters are:\n" << endmsg;
      theOmegaPiLh->printFitParams(std::cout, theParamVal);
      
      //now dump this intermediate fit result
      
      std::ostringstream theParamFilePathGen; 
      theParamFilePathGen << "./" << theAppParams.getName() << "GenevaLastFitParamOmegaPi0Fit_Lmax" << theAppParams.getLMax() << "_mom" << theAppParams.getPbarMom() << ".txt";
      ofstream outfileGen (theParamFilePathGen.str().c_str());
      theOmegaPiLh->dumpCurrentResult(outfileGen, theParamVal);
      outfileGen.close();

      //now fit with minuit

      bExecFinish = Minuit(theAppParams, theParamVal, theOmegaPiLh, false);

      //now dump the final fit result
      std::ostringstream theParamFilePathMin;
      theParamFilePathMin << "./" << theAppParams.getName() << "LastFitParamOmegaPi0Fit_Mmax" << theAppParams.getLMax() << "_mom" << theAppParams.getPbarMom() << ".txt";
      ofstream outfileMin (theParamFilePathMin.str().c_str());
      theOmegaPiLh->dumpCurrentResult(outfileMin, theParamVal);
      outfileMin.close();

      if(bl::indeterminate(bExecFinish)) return 0; // Simply terminate -- this is a client in networked mode that has done its job
      break;
    }

  case ApplicationParameterLS::SpinDensity:
    {
      bExecFinish = calcSpinDensity(theAppParams, theOmegaPiLh, theParamVal);
      break;
    }

  case ApplicationParameterLS::QAmode:
    {
      bExecFinish = QAmode(theAppParams, theParamVal, theOmegaPiLh);
      break;
    }
  default:
    {
      cerr << "Error unknown execution mode selected!" << endl;
    }
  }

  //---------------------------------------------------------------------------
  // Let the audience know
  if (bExecFinish) {
    cout << endl;
    Info << "Fit results:\n" << endmsg;
    Info << "Final logLH=" << theOmegaPiLh->calcLogLh(theParamVal) << "\n" << endmsg;
    Info << "Final fit parameters:\n" << endmsg;
    //       printFitParameters(pbarpToOmegaPi0StatesPtr, finalFitParm);
    theOmegaPiLh->printFitParams(std::cout, theParamVal);
    std::ostringstream theRootFilePath;
    theRootFilePath << "./" << theAppParams.getName() << "OmegaPi0Fit_Lmax" << theAppParams.getLMax() << "_mom" << theAppParams.getPbarMom() << ".root";
    OmegaPiHistLS theHistogrammer(theOmegaPiLh,theParamVal,theRootFilePath.str());
  }

}

