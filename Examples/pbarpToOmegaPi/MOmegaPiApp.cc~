#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiHist.hh"

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

//   std::string piomegaDatFile(theSourcePath+"/Examples/pbarpToOmegaPi/data/510_0600.dat"); 
  std::string piomegaDatFile(theSourcePath+"/Examples/pbarpToOmegaPi/data/510_1940.dat"); 

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

//   std::string piomegaMcFile(theSourcePath+"/Examples/pbarpToOmegaPi/data/mc510_0600.dat"); 
  std::string piomegaMcFile(theSourcePath+"/Examples/pbarpToOmegaPi/data/mc510_1940.dat"); 

  fileNamesMc.push_back(piomegaMcFile);
  CBElsaReader eventReaderMc(fileNamesMc, 3, 0); 
  EventList piOmegaEventsMc;
  eventReaderMc.fillAll(piOmegaEventsMc);
  piOmegaEventsMc.rewind();

  boost::shared_ptr<const OmegaPiEventList> theOmegaPiEventPtr(new OmegaPiEventList(piOmegaEventsData, piOmegaEventsMc));


  OmegaPiHist theHistogrammer(theOmegaPiEventPtr);

//   OmegaPiEventList 
//  PwaEnv::instance().setup(setupFile);

//   ParticleTable* pTable = PwaEnv::instance().particleTable();
//   if (0 == pTable)
//     ErrMsg(fatal) << "getting ParticleTable failed" << endmsg;
//   pTable->print(std::cout);

//   EventList* eventList = PwaEnv::instance().beamEventList();
//   if (0 == eventList)
//     ErrMsg(fatal) << "getting beam EventList failed" << endmsg;

//   ErrMsg(routine) << "Input file has " << eventList->size() << " events. Each event has "
//                   <<  eventList->nextEvent()->size() << " final state particles.\n" << endmsg;
//   eventList->rewind();

//   EventList* mcEventList = PwaEnv::instance().mcEventList();
//   if (0 == mcEventList)
//     ErrMsg(fatal) << "getting MC EventList failed" << endmsg;

//   ErrMsg(routine) << "MC Input file has " << mcEventList->size() << " events. Each event has "
//                   <<  mcEventList->nextEvent()->size() << " final state particles.\n" << endmsg;
//   mcEventList->rewind();

//   Event* anEvent;
//   int evtCount = 0;
//   ErrMsg(routine) << "======== beam events ========" << endmsg;
//   while ((anEvent = eventList->nextEvent()) != 0 && evtCount < 20) {
//     ErrMsg(routine) << "\n" 
//                     << *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
//                     << *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
//                     << *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
//                     << endmsg;
//     ++evtCount;
//   }

//   evtCount = 0;
//   ErrMsg(routine) << "======== MC events ========" << endmsg;
//   while ((anEvent = mcEventList->nextEvent()) != 0 && evtCount < 20) {
//     ErrMsg(routine) << "\n" 
//                     << *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
//                     << *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
//                     << *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
//                     << endmsg;
//     ++evtCount;
//   }



//  std::stringstream dataSpinStrStr(dataSpinStr);
//  int dataSpin=2;
//  dataSpinStrStr >> dataSpin ;

//  ErrMsg(routine) << "dataSpin: " << dataSpin << endmsg;

//  boost::shared_ptr<const EtacToapi0EventList> theEvtListPtr(new EtacToapi0EventList(dataSpin));
 
//  boost::shared_ptr<EtacToapi0Lh> theEtacToapi0LhPtr(new EtacToapi0Lh(theEvtListPtr));
//  MEtacToapi0Fcn fcn(theEtacToapi0LhPtr);

//  MnUserParameters upar;
//  if (dataSpin==2)
//     {  
//       upar.Add("InterMass", 1.6, .1, 2.3, 0.7);
//       upar.Add("InterWidth", 0.04, .01, 0.8, 0.01);
//       upar.Add("spin0", 0.3, .1, 1., 0.);
//       upar.Add("spin1", 0.3, .1, 1., 0.);
//       upar.Add("spin2", 0.3, .1, 1., 0.);
//     }
//  else if (dataSpin==0)
//     {
//      upar.Add("InterMass", 1.1, .1, 2.3, 0.6);
//      upar.Add("InterWidth", 0.04, .01, 0.8, 0.01);
//      upar.Add("spin0", 0.3, .1, 1., 0.);
//      upar.Add("spin1", 0.3, .1, 1., 0.);
//      upar.Add("spin2", 0.3, .1, 1., 0.);
//     }
//  else 
//    {
//      ErrMsg(fatal) << "initialization of the MnUserParameters failed" << endmsg;
//    }

//  MnMigrad migrad(fcn, upar);
//  ErrMsg(routine) <<"start migrad "<< endmsg;
//  FunctionMinimum min = migrad();

//  if(!min.IsValid()) {
//    //try with higher strategy
//    ErrMsg(routine) <<"FM is invalid, try with strategy = 2."<< endmsg;
//    MnMigrad migrad2(fcn, min.UserState(), MnStrategy(2));
//    min = migrad2();
//  }

//  ErrMsg(routine) << "migrad.Fval(): " << min.Fval() << endmsg;
 
//  fitParamVal theFitResult;
//  theFitResult.aMass=min.UserState().Value("InterMass");
//  theFitResult.aWidth=min.UserState().Value("InterWidth");
//  theFitResult.cont0spin=min.UserState().Value("spin0");
//  theFitResult.cont1spin=min.UserState().Value("spin1");
//  theFitResult.cont2spin=min.UserState().Value("spin2");

//  EtacToapi0Hist theHistogrammer(theEvtListPtr,theFitResult);

 if (0!=myLogger) delete myLogger;
 return 0;
}

