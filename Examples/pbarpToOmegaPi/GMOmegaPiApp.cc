/**
 * @file GStartProject.cpp
 */

/* Copyright (C) Dr. Ruediger Berlich and Karlsruhe Institute of Technology
 * (University of the State of Baden-Wuerttemberg and National Laboratory
 * of the Helmholtz Association)
 *
 * Contact: info [at] gemfony (dot) com
 *
 * This file is part of the Geneva library, Gemfony scientific's optimization
 * library.
 *
 * Geneva is free software: you can redistribute it and/or modify
 * it under the terms of version 3 of the GNU Affero General Public License
 * as published by the Free Software Foundation.
 *
 * Geneva is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with the Geneva library.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For further information on Gemfony scientific and Geneva, visit
 * http://www.gemfony.com .
 */

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>

// Boost header files go here
#include <boost/lexical_cast.hpp>

// GenEvA header files go here
#include "GRandom.hpp"
#include "GEvolutionaryAlgorithm.hpp"
#include "GMultiThreadedEA.hpp"
#include "GBrokerEA.hpp"
#include "GIndividualBroker.hpp"
#include "GAsioTCPConsumer.hpp"
#include "GAsioTCPClient.hpp"
#include "GAsioHelperFunctions.hpp"

#include "PwaUtils/pbarpStates.hh"
#include "Examples/pbarpToOmegaPi/AbsOmegaPiLh.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiLhGamma.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiLhOmega.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"
// The individual that should be optimized
#include "Examples/pbarpToOmegaPi/GOmegaPiIndividual.hh"

// Declares a function to parse the command line
#include "Examples/pbarpToOmegaPi/GMArgumentParser.hh"

// Information retrieval and printing
#include "Examples/pbarpToOmegaPi/GInfoFunction.hh"

#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiHist.hh"

#include "Examples/pbarpToOmegaPi/OmegaPiData.hh"

#include "Setup/PwaEnv.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "Event/CBElsaReader.hh"
#include "Particle/PdtParser.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Examples/pbarpToOmegaPi/GMArgumentParser.hh"

#include "Examples/pbarpToOmegaPi/MOmegaPiFcn.hh"

#include "qft++/topincludes/tensor.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

#include "boost/date_time/posix_time/posix_time.hpp"

#include <complex>

#include "Examples/pbarpToOmegaPi/spindensityfitparameter.hh"

#include "Examples/pbarpToOmegaPi/spindensityhist.hh"

using namespace Gem::GenEvA;
using namespace Gem::Util;
using namespace ROOT::Minuit2;
using namespace std;
using namespace boost::posix_time;

inline void printFitParameters(boost::shared_ptr<const pbarpToOmegaPi0States> pbarpToOmegaPi0StatesPtr,
                               OmegaPiData::fitParamVal &theParamVal)
{
          //  print fit paramss
        std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaSinglet=pbarpToOmegaPi0StatesPtr->jpclsSinglet();
        std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0=pbarpToOmegaPi0StatesPtr->jpclsTriplet0();
        std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1=pbarpToOmegaPi0StatesPtr->jpclsTriplet1();
        std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

        Info << "***fit parameter singlet states*** " <<endmsg;  
        for ( itJPCLS=JPCLSOmegaSinglet.begin(); itJPCLS!=JPCLSOmegaSinglet.end(); ++itJPCLS){
          Info << (*itJPCLS)->name()<< endmsg;
          std::pair<double, double> tmpParam=theParamVal.omegaProdSinglet[(*itJPCLS)];
          Info <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
        }
        Info << "***fit parameter triplet m=0 states*** " <<endmsg;  
        for ( itJPCLS=JPCLSOmegaTriplet0.begin(); itJPCLS!=JPCLSOmegaTriplet0.end(); ++itJPCLS){
          Info << (*itJPCLS)->name()<< endmsg;
          std::pair<double, double> tmpParam=theParamVal.omegaProdTriplet0[(*itJPCLS)];
          Info <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
        }
        Info << "***fit parameter triplet m=1 states*** " <<endmsg;  
        for ( itJPCLS=JPCLSOmegaTriplet1.begin(); itJPCLS!=JPCLSOmegaTriplet1.end(); ++itJPCLS){
          Info << (*itJPCLS)->name()<< endmsg;
          std::pair<double, double> tmpParam=theParamVal.omegaProdTriplet1[(*itJPCLS)];
          Info <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
        }
        Info << endmsg;  
}

//This function constructs the path to the file.
inline void constructPath(const string &thePrefix, const unsigned pbarMom, string &outFilePath)
{
  stringstream sstrDatFile; //String Stream for die construction of the path to the parameter File;
  sstrDatFile << thePrefix; 
  sstrDatFile.width(4);
  sstrDatFile.fill('0');
  sstrDatFile << right << pbarMom << ".dat";
  outFilePath = sstrDatFile.str();
}

//This function checks if the file in the path theFilePath exists
inline bool checkFileExist(const string &theFilePath)
{
  ifstream datChk(theFilePath.c_str());
  if (datChk) { return true; } 
  else { return false; }
}

inline const string PrintJPLCS(std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmS,
                       const string &theSuffix
                      )
{   
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;

  ostringstream theOutStream;    
  for ( it=fitParmS.begin(); it!=fitParmS.end(); ++it)
  {
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;
    string strName = theJPCLS->name()+theSuffix;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    theOutStream << strName << " " << theMag << " " << thePhi << endl;
  }
  return theOutStream.str(); 
}

inline const string  PrintFinalFitParam(OmegaPiData::fitParamVal &finalFitParm)
{

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmSinglet=finalFitParm.omegaProdSinglet;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet0=finalFitParm.omegaProdTriplet0;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet1=finalFitParm.omegaProdTriplet1;
  
  ostringstream theOutStream;
  
  theOutStream << PrintJPLCS(fitParmSinglet,"S")
               << PrintJPLCS(fitParmTriplet0,"T0")
               << PrintJPLCS(fitParmTriplet1,"T1");
  
  return theOutStream.str();
}

inline bool GenEvA(boost::shared_ptr<const OmegaPiEventList> &theOmegaPiEventPtr,
		   boost::shared_ptr<const pbarpToOmegaPi0States> &pbarpToOmegaPi0StatesPtr,
		   ApplicationParameter &theAppParams,
		   OmegaPiData::fitParamVal &finalFitParm,
		   boost::shared_ptr<AbsOmegaPiLh> &finalOmegaPiLh
		  )
{
   Info << "GenEvA fit start.\n" << endmsg;
  // Create the first set of parent individuals. Initialization of parameters is done randomly.
  std::vector<boost::shared_ptr<GOmegaPiIndividual> > parentIndividuals;
  for(std::size_t p = 0 ; p<theAppParams.getNParents(); p++) {
    boost::shared_ptr<GOmegaPiIndividual> gdii_ptr( new GOmegaPiIndividual(finalOmegaPiLh) );
    gdii_ptr->setProcessingCycles(theAppParams.getProcessingCycles());

    parentIndividuals.push_back(gdii_ptr);
  }

  //***************************************************************************
  // Create an instance of our optimization monitor, telling it to output information in given intervals
  std::ostringstream theGenResultPath;
  theGenResultPath << "./" << theAppParams.getName() << "result" << "OmegaPi0Fit_jmax" << theAppParams.getJMax() << "_mom" << theAppParams.getPbarMom() << ".C";

  std::ofstream resultSummary(theGenResultPath.str().c_str()); //"./result.C");
  boost::shared_ptr<optimizationMonitor> om(new optimizationMonitor(theAppParams.getNParents(), resultSummary));

  //***************************************************************************
  // We can now start creating populations. We refer to them through the base class

  // This smart pointer will hold the different population types
  boost::shared_ptr<GEvolutionaryAlgorithm> pop_ptr;

  // Create the actual populations
  switch (theAppParams.getParallelizationMode()) {
    //-----------------------------------------------------------------------------------------------------
  case 0: // Serial execution
    // Create an empty population
    pop_ptr = boost::shared_ptr<GEvolutionaryAlgorithm>(new GEvolutionaryAlgorithm());
    break;

    //-----------------------------------------------------------------------------------------------------
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

    //-----------------------------------------------------------------------------------------------------
  case 2: // Networked execution (server-side)
    {
      // Create a network consumer and enrol it with the broker
      boost::shared_ptr<GAsioTCPConsumer> gatc(new GAsioTCPConsumer(theAppParams.getPort()));
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


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Now we have suitable populations and can fill them with data

  // Add individuals to the population
  for(std::size_t p = 0 ; p<theAppParams.getNParents(); p++) {
    pop_ptr->push_back(parentIndividuals[p]);
  }

  // Specify some general population settings
  pop_ptr->setPopulationSize(theAppParams.getPopulationSize(),theAppParams.getNParents());
  pop_ptr->setMaxIteration(theAppParams.getMaxIterations());
  pop_ptr->setMaxTime(boost::posix_time::minutes(theAppParams.getMaxMinutes()));
  pop_ptr->setReportIteration(theAppParams.getReportIteration());
  pop_ptr->setRecombinationMethod(theAppParams.getRScheme());
  pop_ptr->setSortingScheme(theAppParams.getSmode());
  pop_ptr->registerInfoFunction(boost::bind(&optimizationMonitor::informationFunction, om, _1, _2));
  
  // Do the actual optimization
  pop_ptr->optimize();

  //--------------------------------------------------------------------------------------------

  // Make sure we close the result file
  resultSummary.close();

  boost::shared_ptr<GOmegaPiIndividual> bestIndividual_ptr=pop_ptr->getBestIndividual<GOmegaPiIndividual>();
   assert(bestIndividual_ptr->getFitParams(finalFitParm));
  finalOmegaPiLh=bestIndividual_ptr->getOmegaPiLhPtr();

  Info << "GenEvA done.\n" << endmsg;

  return true;


}


inline bool Minuit(boost::shared_ptr<const OmegaPiEventList> &theOmegaPiEventPtr,
		             boost::shared_ptr<const pbarpToOmegaPi0States> &pbarpToOmegaPi0StatesPtr,
		             ApplicationParameter &theAppParams,
		             OmegaPiData::fitParamVal &finalFitParm,
		             boost::shared_ptr<AbsOmegaPiLh> &finalOmegaPiLh
		            )
{
   Info << "Minuit fit start.\n" << endmsg;
    

    // get pbarpToOmegaPi0States pointer back
    boost::shared_ptr<const pbarpToOmegaPi0States> theOmegaPi0StatesPtr=finalOmegaPiLh->omegaPi0States();  
    
    theOmegaPi0StatesPtr->print(std::cout);


    MOmegaPiFcn mOmegaPiFcn(finalOmegaPiLh);

   MnUserParameters upar;
   minuitStartParam theUserPar;

   if(theAppParams.getAppExecMode() == ApplicationParameter::Minuit)
   {
      if (!theAppParams.getPathStartParamFile().empty())
      {
         ifstream theFile(theAppParams.getPathStartParamFile().c_str());   
         if (theFile)
         {
            Info << "Using start parameters from file " << theAppParams.getPathStartParamFile() << "\n" << endmsg;
            theUserPar.ParseStream(theFile);
            mOmegaPiFcn.setMnUsrParams(upar,theUserPar);
         }
         else 
         {
            Info << "Start parameter file " << theAppParams.getPathStartParamFile() << " not found!" << endmsg;
            return false;
         }
      }
      else mOmegaPiFcn.setMnUsrParams(upar);
   }
   else if(theAppParams.getAppExecMode() == ApplicationParameter::GenToMinuit) mOmegaPiFcn.setMnUsrParams(upar,finalFitParm);
   
  
   MnMigrad migrad(mOmegaPiFcn, upar);
   Info <<"start migrad "<< endmsg;
   FunctionMinimum min = migrad();

  if(!min.IsValid()) {
    //try with higher strategy
    Info <<"FM is invalid, try with strategy = 2."<< endmsg;
    MnMigrad migrad2(mOmegaPiFcn, min.UserState(), MnStrategy(2));
    min = migrad2();
  }

  MnUserParameters finalUsrParameters=min.UserParameters();
  const std::vector<double> finalParamVec=finalUsrParameters.Params();

  cout << endl << "Errors:" << endl;
  const vector<double> finalParamErrors=finalUsrParameters.Errors();
  
  vector<double>::const_iterator it;
  int i;
  for (it = finalParamErrors.begin(), i=0; it < finalParamErrors.end(); it++,i++) 
  {
    string strName =finalUsrParameters.GetName(i);
    cout << "Name:" << strName << " ";
    cout << "Index:" << finalUsrParameters.Index(strName) << " ";
    cout << *it << endl;
  }
  
  mOmegaPiFcn.setFitParamVal(finalFitParm, finalParamVec);
  
  Info << "Minuit done.\n" << endmsg;

  return true;
}

inline bool QAmode(boost::shared_ptr<const OmegaPiEventList> &theOmegaPiEventPtr,
		             boost::shared_ptr<const pbarpToOmegaPi0States> &pbarpToOmegaPi0StatesPtr,
		             ApplicationParameter &theAppParams,
		             OmegaPiData::fitParamVal &finalFitParm,
		             boost::shared_ptr<AbsOmegaPiLh> &finalOmegaPiLh
		            )
{
   Info << "QA mode start.\n" << endmsg;
    

    // get pbarpToOmegaPi0States pointer back
    boost::shared_ptr<const pbarpToOmegaPi0States> theOmegaPi0StatesPtr=finalOmegaPiLh->omegaPi0States();  
    
    theOmegaPi0StatesPtr->print(std::cout);


    MOmegaPiFcn mOmegaPiFcn(finalOmegaPiLh);

    MnUserParameters upar;
    minuitStartParam theUserPar;
    

    if (!theAppParams.getPathStartParamFile().empty())
      {
	ifstream theFile(theAppParams.getPathStartParamFile().c_str());   
	if (theFile)
	  {
	    Info << "Using start parameters from file " << theAppParams.getPathStartParamFile() << "\n" << endmsg;
	    theUserPar.ParseStream(theFile);
	    mOmegaPiFcn.setMnUsrParams(upar,theUserPar);
	  }
	else 
	  {
	    Info << "Start parameter file " << theAppParams.getPathStartParamFile() << " not found!" << endmsg;
	    return false;
	  }
      }
	else 
	  {
	    Info << "Start parameter file has to be set!" << endmsg;
	    return false;
	  }
      
	mOmegaPiFcn.setFitParamVal(finalFitParm, upar.Params());
	
	Info << "QA mode done.\n" << endmsg;
	
	return true;
      }

inline bool calcSpinDensity(ApplicationParameter &theAppParams, 
                            boost::shared_ptr<AbsOmegaPiLh> finalOmegaPiLh, 
                            boost::shared_ptr<const pbarpToOmegaPi0States> pbarpToOmegaPi0StatesPtr,
                            boost::shared_ptr<const OmegaPiEventList> theOmegaPiEventPtr)
{
  Info << "Starting spin density calculation." << "\n" << endmsg;
  SpinDensityFitParameter theUserPar;
  OmegaPiData::fitParamVal theParamVal;
        
  if (!theAppParams.getPathStartParamFile().empty())
  {
    ifstream theFile(theAppParams.getPathStartParamFile().c_str());   
    if (theFile)
    {
      Info << "Using start parameters from file " << theAppParams.getPathStartParamFile() << "\n" << endmsg;
      theUserPar.ParseStream(theFile);
      theUserPar.getFitParamVal(theParamVal,finalOmegaPiLh->omegaPi0States());
    }
    else 
    {
      cerr << "Fehler" << endl;
      Info << "Start parameter file " << theAppParams.getPathStartParamFile() << " not found!" << endmsg;
      return false;
    }
  }
  else 
  {
    cerr << "Start parameter file not specified!" << endl;
    return false;
  }
        
  Info << "Using following fit parameter:\n" << endmsg;
  printFitParameters(pbarpToOmegaPi0StatesPtr, theParamVal);
  std::ostringstream theSpinDensityRootFile;
        
  if (theAppParams.getCalcAllSpindensity())
  {
    Info << "Calculating all spin density elements.\n" << endmsg;
    theSpinDensityRootFile << "./" << theAppParams.getName() << "SpinDensity" << "_jmax" << theAppParams.getJMax() << "_mom" << theAppParams.getPbarMom() << ".root";

    SpinDensityHist theSpinDensityHist(theSpinDensityRootFile.str(),theOmegaPiEventPtr->getMcVecs(),theParamVal);
    theSpinDensityHist.createHistograms();
  }
  else
  {
    Info << "Calculating spin density elements for M=" << theAppParams.getM() << " M'=" << theAppParams.getM_() << ".\n" << endmsg;
    theSpinDensityRootFile << "./" << theAppParams.getName() << "SpinDensity_M1" << theAppParams.getM() << "_M2" << theAppParams.getM_() << "_jmax" << theAppParams.getJMax() << "_mom" << theAppParams.getPbarMom() << ".root";
       
    SpinDensityHist theSpinDensityHist(theSpinDensityRootFile.str(),theOmegaPiEventPtr->getMcVecs(),theParamVal);
    theSpinDensityHist.createHistogram(theAppParams.getM(),theAppParams.getM_(), true);
    theSpinDensityHist.createHistogram(theAppParams.getM(),theAppParams.getM_(), false);
  }
  Info << "Spin density calculation done." << "\n" << endmsg;
  return true;
}

inline void removeEvents(EventList &piOmegaEventsData, int nEventsToRemove, bool bRemoveFromEnd)
{
  if (nEventsToRemove != 0)
  {
    stringstream strMsg;
    strMsg << "Removing " << nEventsToRemove << " events";
    if (bRemoveFromEnd)
    {
      strMsg << " from the end of the event list."  << endl;
      piOmegaEventsData.removeEvents(piOmegaEventsData.size()-nEventsToRemove,piOmegaEventsData.size());
    }
    else
    {
      strMsg << " from the beginning of the event list." << endl;
      piOmegaEventsData.removeEvents(0,piOmegaEventsData.size()-nEventsToRemove);
    }
    Info << strMsg.str() << endmsg;
    Info << "\nEvent list now contains " << piOmegaEventsData.size() << " events. Each event has "
        <<  piOmegaEventsData.nextEvent()->size() << " final state particles.\n" << endmsg;
  }
}

/************************************************************************************************/
/**
 * The main function.
 */
int main(int argc, char **argv)
{
  try
  {

    cout << "Combined GenEvA and Minuit fit program " << argv[0] << endl << endl;

    static ApplicationParameter theAppParams(argc, argv);
    
    ptime startTime = second_clock::local_time();
    cout << "Start time : " << startTime << endl;
    
    // Random numbers are our most valuable good. Set the number of threads
    GRANDOMFACTORY->setNProducerThreads(theAppParams.getNProducerThreads());
    GRANDOMFACTORY->setArraySize(theAppParams.getArraySize());

  switch(theAppParams.getErrLogMode()) {
	case ApplicationParameter::debug :
	ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
	break;
      case ApplicationParameter::trace :
	ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
	break;
      case ApplicationParameter::routine :
	ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
	break;
      case ApplicationParameter::warning :
	ErrLogger::instance()->setLevel(log4cpp::Priority::WARN);
	break;
      case ApplicationParameter::error :
	ErrLogger::instance()->setLevel(log4cpp::Priority::ERROR);
	break;
      case ApplicationParameter::alert :
	ErrLogger::instance()->setLevel(log4cpp::Priority::ALERT);
	break;
      default: 
	ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
      }

      switch(theAppParams.getAppExecMode())
      {
	     case ApplicationParameter::GenEvA:
	        Info << "Using minimization algorithm: " << "GenEvA" << "\n" << endmsg;
	        break;
	     case ApplicationParameter::Minuit:
          Info << "Using minimization algorithm: " << "Minuit" << "\n" << endmsg;
	        break;
        case ApplicationParameter::GenToMinuit:
          Info << "Using minimization algorithm: " << "First GenEvA then Minuit with final fit parameters from GenEvA" << "\n" << endmsg;
           break;
        case ApplicationParameter::SpinDensity:
          Info << "Calculating Spin Density" << "\n" << endmsg;
          break;
        case ApplicationParameter::QAmode:
          Info << "using QA mode" << "\n" << endmsg;
          break;

      }

      string theCfgFile = theAppParams.getConfigFile();
      Info << "The path to config file is " << theCfgFile << "\n" << endmsg;
      //***************************************************************************
      // If this is a client in networked mode, we can just start the listener and
      // return when it has finished
    if(theAppParams.getParallelizationMode()==2 && !theAppParams.getServerMode()) {
	boost::shared_ptr<GAsioTCPClient> p(new GAsioTCPClient(theAppParams.getIp(), boost::lexical_cast<std::string>(theAppParams.getPort())));

	p->setMaxStalls(0); // An infinite number of stalled data retrievals
	p->setMaxConnectionAttempts(100); // Up to 100 failed connection attempts

	// Prevent return of unsuccessful mutation attempts to the server
	p->returnResultIfUnsuccessful(theAppParams.getReturnRegardless());

	// Start the actual processing loop
	p->run();

	return 0;
      }
      //***************************************************************************

      Info << "Maximum spin content: " << theAppParams.getJMax() << endmsg;
      Info << "pbar momentum: " << theAppParams.getPbarMom()   << endmsg;
      
  std::string piomegaDatFile;
  std::string piomegaMcFile; 
  
  constructPath(theAppParams.getSourcePath()+"/Examples/pbarpToOmegaPi/data/510_",theAppParams.getPbarMom(),piomegaDatFile);
  if (checkFileExist(piomegaDatFile))
  {
    Info << "Using Data file " << piomegaDatFile << endmsg;
  }
  else
  {
    Alert <<"Data file for pbarMom= " << theAppParams.getPbarMom() << " not available !"  << endmsg;
    Alert << "File " << piomegaDatFile << " is missing !" << endmsg;
    exit(1);
  }

  constructPath(theAppParams.getSourcePath()+ "/Examples/pbarpToOmegaPi/data/mc510_",theAppParams.getPbarMom(),piomegaMcFile);
  if (checkFileExist(piomegaMcFile))
  {
    Info << "Using Monte Carlo file " << piomegaMcFile << endmsg;
  }
  else
  {
    Alert <<"Monte Carlo file for pbarMom= " << theAppParams.getPbarMom() << "not available !" << endmsg;
    Alert << "File " << piomegaMcFile << " is missing !" << endmsg;
    exit(1);
  }
  
  Info << "data file: " << piomegaDatFile << endmsg;
  Info << "mc file: " << piomegaMcFile << endmsg;

      ParticleTable pTable;
      PdtParser parser;
      std::string pdtFile(theAppParams.getSourcePath()+"/Particle/pdt.table");
      if (!parser.parse(pdtFile, pTable)) {
	Alert << "Error: could not parse " << pdtFile << endmsg;
	exit(1);
      }

      std::vector<std::string> fileNames;

      fileNames.push_back(piomegaDatFile);
      CBElsaReader eventReader(fileNames, 3, 0); 
      EventList piOmegaEventsData;
      eventReader.fillAll(piOmegaEventsData);

      Info << "\nFile has " << piOmegaEventsData.size() << " events. Each event has "
          <<  piOmegaEventsData.nextEvent()->size() << " final state particles.\n" << endmsg;

      if(theAppParams.getNumEventsRemove() > 0) 
      {
        removeEvents(piOmegaEventsData,theAppParams.getNumEventsRemove(),theAppParams.getRemoveEventsFromEnd());
      }
      else if(theAppParams.getNumEventsRed() > 0)
      {
        int nEventsRed = theAppParams.getNumEventsRed();
        int nEventListSize = piOmegaEventsData.size();
        
        Info << "Reducing number of the events in the event list to " 
            << nEventsRed << " events." << endmsg;
        
        int nNewSize = nEventListSize-nEventsRed;
        
        if (nNewSize > 0 && nNewSize != nEventListSize)
        {
          removeEvents(piOmegaEventsData, nNewSize, theAppParams.getRemoveEventsFromEnd());
        }
        else
        {
          Alert << "Could not reduce number of the parameters because number to which they should be"
                  << " reduced " << "(" << nEventsRed << ")" 
                  <<" to is to great or same as the actual number of events in the list.\n" 
                  << endmsg;
          exit(-1);
        }
      }
      
      if (!piOmegaEventsData.findParticleTypes(pTable))	Warning << "could not find all particles" << endmsg;

      piOmegaEventsData.rewind();


      Event* anEvent;
      int evtCount = 0;
      while ((anEvent = piOmegaEventsData.nextEvent()) != 0 && evtCount < 20) {
	Info << "\n" 
			<< *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
			<< *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
			<< *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
			<< endmsg;
	++evtCount;
      }
      piOmegaEventsData.rewind();

      std::vector<std::string> fileNamesMc;

      fileNamesMc.push_back(piomegaMcFile);
      CBElsaReader eventReaderMc(fileNamesMc, 3, 0); 
      EventList piOmegaEventsMc;
      eventReaderMc.fillAll(piOmegaEventsMc);
      piOmegaEventsMc.rewind();

      boost::shared_ptr<const OmegaPiEventList> theOmegaPiEventPtr(new OmegaPiEventList(piOmegaEventsData, piOmegaEventsMc, theAppParams.getJMax(),  theAppParams.getPbarMom()));

      boost::shared_ptr<pbarpStates> pbarpStatesPtr(new pbarpStates(theAppParams.getJMax()));

      boost::shared_ptr<const pbarpToOmegaPi0States> pbarpToOmegaPi0StatesPtr(new pbarpToOmegaPi0States(pbarpStatesPtr));
    

      OmegaPiData::fitParamVal finalFitParm;


      boost::shared_ptr<AbsOmegaPiLh> finalOmegaPiLh= boost::shared_ptr<AbsOmegaPiLh>(new OmegaPiLhGamma(theOmegaPiEventPtr, pbarpToOmegaPi0StatesPtr));

      ApplicationParameter::enExecMode theExecMode = theAppParams.getAppExecMode();
      bool bExecFinish=false;      
      if(theExecMode == ApplicationParameter::GenEvA)
      {
	       bExecFinish = GenEvA(theOmegaPiEventPtr,pbarpToOmegaPi0StatesPtr,theAppParams,finalFitParm,finalOmegaPiLh);
          
          std::ostringstream theParamFilePath;
          theParamFilePath << "./" << theAppParams.getName() << "LastFitParamOmegaPi0Fit_jmax" << theAppParams.getJMax() << "_mom" << theAppParams.getPbarMom() << ".txt";
          ofstream outfile (theParamFilePath.str().c_str());
          outfile << PrintFinalFitParam(finalFitParm);
          outfile.close();
      }
      else if(theExecMode == ApplicationParameter::Minuit)
      {
          bExecFinish = Minuit(theOmegaPiEventPtr,pbarpToOmegaPi0StatesPtr,theAppParams,finalFitParm,finalOmegaPiLh);
      }
      else if(theExecMode == ApplicationParameter::GenToMinuit)
      {
            bExecFinish = GenEvA(theOmegaPiEventPtr,pbarpToOmegaPi0StatesPtr,theAppParams,finalFitParm,finalOmegaPiLh);
            bExecFinish = Minuit(theOmegaPiEventPtr,pbarpToOmegaPi0StatesPtr,theAppParams,finalFitParm,finalOmegaPiLh);
      }
      else if(theExecMode == ApplicationParameter::SpinDensity)
      {
        if (!calcSpinDensity(theAppParams, finalOmegaPiLh, pbarpToOmegaPi0StatesPtr,theOmegaPiEventPtr)) exit(1);
      }
      else if(theExecMode == ApplicationParameter::QAmode)
	{
	   bExecFinish = QAmode(theOmegaPiEventPtr,pbarpToOmegaPi0StatesPtr,theAppParams,finalFitParm,finalOmegaPiLh);
	}
      else
      {
         cerr << "Error unknown execution mode selected!" << endl;
      }

      if (bExecFinish) 
      {
        cout << endl;
        Info << "Fit results:\n" << endmsg;
        Info << "Final logLH=" << finalOmegaPiLh->calcLogLh(finalFitParm) << "\n" << endmsg;
        Info << "Final fit parameters:\n" << endmsg;
        printFitParameters(pbarpToOmegaPi0StatesPtr, finalFitParm);

        std::ostringstream theRootFilePath;
        theRootFilePath << "./" << theAppParams.getName() << "OmegaPi0Fit_jmax" << theAppParams.getJMax() << "_mom" << theAppParams.getPbarMom() << ".root";
        OmegaPiHist theHistogrammer(finalOmegaPiLh,finalFitParm,theRootFilePath.str());
      }

      ptime endTime = second_clock::local_time();

      cout << "All done." << endl;
      cout << "End time : " << endTime << endl;
      cout << "Duration : " << (endTime - startTime) << endl;

      return 0;
  }
  catch(...)
  {
    cerr << "Program terminated because of an error!" << endl;
    exit(-1);
  }
}
