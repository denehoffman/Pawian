#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiHist.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiData.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiLh.hh"
#include "Examples/pbarpToOmegaPi/MOmegaPiFcn.hh"

#include "Setup/PwaEnv.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "Event/CBElsaReader.hh"
#include "Particle/PdtParser.hh"
#include "qft++/topincludes/tensor.hh"

#include "ErrLogger/ErrLineLog.hh"
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

    std::cout << "This application is a simple PWA fit for the decay chain\n"
	      << "eta_c -> intermediate + pi0  ; intermediate -> pi0 + eta\n"
	      << "It makes use of the Covariant Tensor Formalism\n"
	      << "The fit determines whether the intermediate resonance is a spin 0,1 or 2 particle\n"
	      << "In addition the mass and width of the intermediate resonance will be fitted with a simple Breit-Wigner\n"
              << "To start the application with data containing  an intermediate resonance with Spin=0 and mass 0.98 GeV, type: ./MEtacToapi0FitApp -d 0\n"
              << "To start the application with data containing  an intermediate resonance with Spin=2 and mass 1.32 GeV, type: ./MEtacToapi0FitApp -d 2\n"
              << "with the flag -msg <errorLogMode> you can choose the mode for the error logger\n"
	      << "i.e. with './MEtacToapi0FitApp -d 2 -msg debugging' you start the fit with Spin=2 in the debugging mode for the error logger\n"  
	      << std::endl;
    return 0;
  }

  
  std::string msgModeStr="default";

  // decode arguments
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-msg"){
      optind++;
      msgModeStr = __argv[optind];
      found=true;
    }
    if (!found){
      ErrMsg(warning) << "Unknown switch: " 
            << __argv[optind] << endmsg;
      optind++;
    }
    
    while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
  }

  ErrLineLog* myLogger=0;
  if(msgModeStr == "debugging") myLogger= new ErrLineLog(ErrLog::debugging);
  else if(msgModeStr == "trace") myLogger= new ErrLineLog(ErrLog::trace);
  else if(msgModeStr == "routine") myLogger= new ErrLineLog(ErrLog::routine);
  else if(msgModeStr == "warning")  myLogger= new ErrLineLog(ErrLog::warning);
  else if(msgModeStr == "error")    myLogger= new ErrLineLog(ErrLog::error); 
  else {
    myLogger= new ErrLineLog(ErrLog::debugging);
    ErrMsg(warning) << "ErrorLogger not (properly) set -> Use mode 'ErrLog::debugging' " << endmsg;  
  }

  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");
  
  ParticleTable pTable;
  PdtParser parser;
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");
  if (!parser.parse(pdtFile, pTable)) {
    ErrMsg(fatal) << "Error: could not parse " << pdtFile << endmsg;
    exit(1);
  }
  

  std::vector<std::string> fileNames;

  std::string piomegaDatFile(theSourcePath+"/Examples/pbarpToOmegaPi/data/510_0600.dat"); 
//   std::string piomegaDatFile(theSourcePath+"/Examples/pbarpToOmegaPi/data/510_1940.dat"); 

  fileNames.push_back(piomegaDatFile);
  CBElsaReader eventReader(fileNames, 3, 0); 
  EventList piOmegaEventsData;
  eventReader.fillAll(piOmegaEventsData);

  if (!piOmegaEventsData.findParticleTypes(pTable))
    ErrMsg(warning) << "could not find all particles" << endmsg;

  ErrMsg(routine) << "\nFile has " << piOmegaEventsData.size() << " events. Each event has "
                  <<  piOmegaEventsData.nextEvent()->size() << " final state particles.\n" << endmsg;
  piOmegaEventsData.rewind();

  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = piOmegaEventsData.nextEvent()) != 0 && evtCount < 20) {
//     double gammaE=anEvent->p4(1)->E()- anEvent->p4(2)->E();
//     double gammaPx=anEvent->p4(1)->Px()- anEvent->p4(2)->Px();
//     double gammaPy=anEvent->p4(1)->Py()- anEvent->p4(2)->Py();
//     double gammaPz=anEvent->p4(1)->Pz()- anEvent->p4(2)->Pz();
//     anEvent->addParticle(gammaE, gammaPx, gammaPy, gammaPz);
    
    ErrMsg(routine) << "\n" 
                    << *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
                    << *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
                    << *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
                    << endmsg;
    ++evtCount;
  }
  piOmegaEventsData.rewind();


  std::vector<std::string> fileNamesMc;

  std::string piomegaMcFile(theSourcePath+"/Examples/pbarpToOmegaPi/data/mc510_0600.dat"); 
//   std::string piomegaMcFile(theSourcePath+"/Examples/pbarpToOmegaPi/data/mc510_1940.dat"); 

  fileNamesMc.push_back(piomegaMcFile);
  CBElsaReader eventReaderMc(fileNamesMc, 3, 0); 
  EventList piOmegaEventsMc;
  eventReaderMc.fillAll(piOmegaEventsMc);
  piOmegaEventsMc.rewind();

  boost::shared_ptr<const OmegaPiEventList> theOmegaPiEventPtr(new OmegaPiEventList(piOmegaEventsData, piOmegaEventsMc, 3));

  boost::shared_ptr<pbarpStates> pbarpStatesPtr(new pbarpStates(3));
  boost::shared_ptr<pbarpToOmegaPi0States> pbarpToOmegaPi0StatesPtr(new pbarpToOmegaPi0States(pbarpStatesPtr));

  boost::shared_ptr<OmegaPiLh> theOmegaLhPtr(new OmegaPiLh(theOmegaPiEventPtr, pbarpToOmegaPi0StatesPtr));

  // get pbarpToOmegaPi0States pointer back
  boost::shared_ptr<const pbarpToOmegaPi0States> theOmegaPi0StatesPtr=theOmegaLhPtr->omegaPi0States();  
  
  theOmegaPi0StatesPtr->print(std::cout);


  MOmegaPiFcn mOmegaPiFcn(theOmegaLhPtr);
  
  MnUserParameters upar;
  mOmegaPiFcn.setMnUsrParams(upar);

  MnMigrad migrad(mOmegaPiFcn, upar);
  ErrMsg(routine) <<"start migrad "<< endmsg;
  FunctionMinimum min = migrad();

 if(!min.IsValid()) {
   //try with higher strategy
   ErrMsg(routine) <<"FM is invalid, try with strategy = 2."<< endmsg;
   MnMigrad migrad2(mOmegaPiFcn, min.UserState(), MnStrategy(2));
   min = migrad2();
 }

 MnUserParameters finalUsrParameters=min.UserParameters();
 const std::vector<double> finalParamVec=finalUsrParameters.Params();

 OmegaPiData::fitParamVal finalFitParams;
 mOmegaPiFcn.setFitParamVal(finalFitParams, finalParamVec);


 OmegaPiHist theHistogrammer(theOmegaLhPtr, finalFitParams);
//   // now fill the fir parameter
//   OmegaPiData::fitParamVal theFitParameter; 
 
//   std::vector< boost::shared_ptr<const JPCLS> > allJPCLSStates=theOmegaPi0StatesPtr->jpclsStates();
//   std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

//   double counter=0.;
//   std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaSinglet=theOmegaPi0StatesPtr->jpclsSinglet();
//   for ( itJPCLS=JPCLSOmegaSinglet.begin(); itJPCLS!=JPCLSOmegaSinglet.end(); ++itJPCLS){
//     //now fill the fitParameterMap
//     std::pair <double,double> tmpParameter=make_pair(counter,0.0);
//     theFitParameter.omegaProdSinglet[(*itJPCLS)]=tmpParameter;
//     counter++; 
//   }

//   std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0=theOmegaPi0StatesPtr->jpclsTriplet0();
//   for ( itJPCLS=JPCLSOmegaTriplet0.begin(); itJPCLS!=JPCLSOmegaTriplet0.end(); ++itJPCLS){
//     //now fill the fitParameterMap
//     std::pair <double,double> tmpParameter=make_pair(counter,0.0);
//     theFitParameter.omegaProdTriplet0[(*itJPCLS)]=tmpParameter;
//     counter++; 
//   }

//   std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1=theOmegaPi0StatesPtr->jpclsTriplet1();
//   for ( itJPCLS=JPCLSOmegaTriplet1.begin(); itJPCLS!=JPCLSOmegaTriplet1.end(); ++itJPCLS){
//     //now fill the fitParameterMap
//     std::pair <double,double> tmpParameter=make_pair(counter,0.0);
//     theFitParameter.omegaProdTriplet1[(*itJPCLS)]=tmpParameter;
//     counter++; 
//   }

//   //now print the fitparameter
//   std::cout << "\n The fit parameter for omega production are (singlet states):" << std::endl;
//   for ( itJPCLS=JPCLSOmegaSinglet.begin(); itJPCLS!=JPCLSOmegaSinglet.end(); ++itJPCLS){
//     std::cout << (*itJPCLS)->name();
//     std::pair<double, double> tmpParam=theFitParameter.omegaProdSinglet[(*itJPCLS)];
//     std::cout <<"\t" << tmpParam.first <<"\t" << tmpParam.second  << std::endl;
//   }

//   std::cout << "\n The fit parameter for omega production are (triplet0 states):" << std::endl;
//   for ( itJPCLS=JPCLSOmegaTriplet0.begin(); itJPCLS!=JPCLSOmegaTriplet0.end(); ++itJPCLS){
//     std::cout << (*itJPCLS)->name();
//     std::pair<double, double> tmpParam=theFitParameter.omegaProdTriplet0[(*itJPCLS)];
//     std::cout <<"\t" << tmpParam.first <<"\t" << tmpParam.second  << std::endl;
//   }

//   std::cout << "\n The fit parameter for omega production are (triplet1 states):" << std::endl;
//   for ( itJPCLS=JPCLSOmegaTriplet1.begin(); itJPCLS!=JPCLSOmegaTriplet1.end(); ++itJPCLS){
//     std::cout << (*itJPCLS)->name();
//     std::pair<double, double> tmpParam=theFitParameter.omegaProdTriplet1[(*itJPCLS)];
//     std::cout <<"\t" << tmpParam.first <<"\t" << tmpParam.second  << std::endl;
//   }

//   OmegaPiHist theHistogrammer(theOmegaPiEventPtr);



 if (0!=myLogger) delete myLogger;
 return 0;
}

