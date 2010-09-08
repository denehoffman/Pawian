#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/MATpbarpToOmegaPi/OmegaPiEventList.hh"
#include "Examples/MATpbarpToOmegaPi/OmegaPiHist.hh"
#include "Examples/MATpbarpToOmegaPi/OmegaPiData.hh"
#include "Examples/MATpbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "Examples/MATpbarpToOmegaPi/OmegaPiLh.hh"
#include "Examples/MATpbarpToOmegaPi/MOmegaPiFcn.hh"

#include "Setup/PwaEnv.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "Event/CBElsaReader.hh"
#include "Particle/PdtParser.hh"
#include "qft++/topincludes/tensor.hh"

#include "ErrLogger/ErrLogger.hh"
//#include "Minuit2/MnUserParameters.h"
//#include "Minuit2/MnMigrad.h"
//#include "Minuit2/FunctionMinimum.h"
//#include "Minuit2/MnMinos.h"
//#include "Minuit2/MnStrategy.h"
#include "Examples/MATpbarpToOmegaPi/PwaIntMinuit.hh"

#include "PwaUtils/pbarpStates.hh"

//#include "Minuit2/MnUserTransformation.h"
using namespace ROOT::Minuit2;

int main(int __argc,char *__argv[]){

  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
		    || strcmp( __argv[1], "--help" ) == 0 ) ){

    Info << "This application is a PWA fit for the reaction\n\n"
	 << "pbar p -> omega + pi0  ; omega -> pi0 + gamma\n\n"
	 << "It makes use of the Helicity Formalism\n"
	 << "The fit determines the content of the initial states formed via the pbar p reaction "
	 << "and the content of the LS combinations for the dacay to omega pi0\n\n"
	 << "For the steering of the application one can use the following flags\n"
	 << "-jmax: maximum spin concidering in the fit (default 3)\n"
	 << "-pbarmom: pbar momentum of the input data (available so far: 600, 1940) (default 600)\n" 
	 << "-msg <errorLogMode> you can choose the mode for the error logger\n\n"
	 << "i.e. with './MpbarpToOmegaPiApp -jmax 5 -pbarmom 1940 -msg debugging' you start the fit with jmax=3, reconstructed data and MCs at 1940 MeV/c in the debugging mode for the error logger\n"  
	 << endmsg;
    return 0;
  }

  
  std::string msgModeStr="default";
  std::string jmaxStr="3";
  std::string pbarmomStr="600";

  // decode arguments
  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-msg"){
      optind++;
      msgModeStr = __argv[optind];
      found=true;
    }
    if (sw=="-jmax"){
      optind++;
      jmaxStr = __argv[optind];
      found=true;
    }
    if (sw=="-pbarmom"){
      optind++;
      pbarmomStr = __argv[optind];
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


  std::stringstream jmaxStrStr(jmaxStr);
  std::stringstream pbarmomStrStr(pbarmomStr);

  unsigned jMax=3;
  jmaxStrStr >> jMax;

  unsigned pbarMom=600;
  pbarmomStrStr >> pbarMom; 
  
  Info << "Maximum spin content: " << jMax << endmsg;
  Info << "pbar momentum: " << pbarMom   << endmsg;

  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");

  std::string piomegaDatFile;
  std::string piomegaMcFile; 
  if(pbarMom==600){
    piomegaDatFile=theSourcePath+"/Examples/pbarpToOmegaPi/data/510_0600.dat";
    piomegaMcFile=theSourcePath+"/Examples/pbarpToOmegaPi/data/mc510_0600.dat";
  }
  else if (pbarMom==1940){
    piomegaDatFile=theSourcePath+"/Examples/pbarpToOmegaPi/data/510_1940.dat";
    piomegaMcFile=theSourcePath+"/Examples/pbarpToOmegaPi/data/mc510_1940.dat";
  }
  else{
    Alert <<"data for pbarMom= " << pbarMom << "not available; use 600 or 1940!!!" << endmsg;
    exit(1);
  }

  Info << "data file: " << piomegaDatFile << endmsg;
  Info << "mc file: " << piomegaMcFile << endmsg;
  
  
  ParticleTable pTable;
  PdtParser parser;
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");
  if (!parser.parse(pdtFile, pTable)) {
    Alert << "Error: could not parse " << pdtFile << endmsg;
    exit(1);
  }
  

  std::vector<std::string> fileNames;

  fileNames.push_back(piomegaDatFile);
  CBElsaReader eventReader(fileNames, 3, 0); 
  EventList piOmegaEventsData;
  eventReader.fillAll(piOmegaEventsData);

  if (!piOmegaEventsData.findParticleTypes(pTable))
    Warning << "could not find all particles" << endmsg;

  Info << "\nFile has " << piOmegaEventsData.size() << " events. Each event has "
                  <<  piOmegaEventsData.nextEvent()->size() << " final state particles.\n" << endmsg;
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

  //boost::shared_ptr<const OmegaPiEventList> theOmegaPiEventPtr(new OmegaPiEventList(piOmegaEventsData, piOmegaEventsMc, jMax,  pbarMom));
  boost::shared_ptr<OmegaPiEventList> theOmegaPiEventPtr = OmegaPiEventList::getList(piOmegaEventsData, piOmegaEventsMc, jMax,  pbarMom);

  boost::shared_ptr<pbarpStates> pbarpStatesPtr(new pbarpStates(jMax));
  //boost::shared_ptr<pbarpToOmegaPi0States> pbarpToOmegaPi0StatesPtr(new pbarpToOmegaPi0States(pbarpStatesPtr));
  boost::shared_ptr<pbarpToOmegaPi0States> pbarpToOmegaPi0StatesPtr = pbarpToOmegaPi0States::getStates(pbarpStatesPtr);

  boost::shared_ptr<OmegaPiLh> theOmegaLhPtr(new OmegaPiLh(pbarpToOmegaPi0StatesPtr));

  // get pbarpToOmegaPi0States pointer back
  boost::shared_ptr<const pbarpToOmegaPi0States> theOmegaPi0StatesPtr=theOmegaLhPtr->omegaPi0States();  
  
  theOmegaPi0StatesPtr->print(std::cout);


  MOmegaPiFcn mOmegaPiFcn(theOmegaLhPtr);
  
  /*MnUserParameters upar;
  mOmegaPiFcn.setMnUsrParams(upar);

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
 const std::vector<double> finalParamVec=finalUsrParameters.Params();*/
 
 PwaIntMinuit test;
 cout << "bla" << endl;
 const std::vector<double> finalParamVec=test.exec(&mOmegaPiFcn);

 OmegaPiData::fitParamVal finalFitParams;
 mOmegaPiFcn.setFitParamVal(finalFitParams, finalParamVec);


 OmegaPiHist theHistogrammer(theOmegaLhPtr, finalFitParams);

//   OmegaPiHist theHistogrammer(theOmegaPiEventPtr);


 return 0;
}

