#include <getopt.h>

#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamHitAndMiss.hh"
#include "Examples/Psi2STo2K2PiGam/AbsPsi2STo2K2PiGamLh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtGenConverter.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamHist.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Event/HepMCEventList.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PspGen/EvtRandomEngine.hh"
#include "PspGen/EvtRandom.hh"
#include "PspGen/EvtSimpleRandomEngine.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_Ascii.h"

Psi2STo2K2PiGamHitAndMiss::Psi2STo2K2PiGamHitAndMiss(boost::shared_ptr<AbsPsi2STo2K2PiGamLh> lhPtr, const param2K2PiGam& theParamVal, boost::shared_ptr<HepMCEventList> hepMCEvtListPtr, std::string histName) :
  _lhPtr(lhPtr)
  ,_psi2STo2K2PiGamHist(new Psi2STo2K2PiGamHist(theParamVal, histName))
{
  
  //first calculate max weight
  double maxWeight=0.;
  hepMCEvtListPtr->rewind();
  HepMC::GenEvent* genEvent;

  while ( (genEvent = hepMCEvtListPtr->nextEvent()) !=0 ){
    Psi2STo2K2PiGamEvtData* dataEvt=Psi2STo2K2PiGamEvtGenConverter::instance()->convertHepMCtoDataEvt(genEvent);

    double currentWeight=lhPtr->calcEvtIntensity(dataEvt, theParamVal);
   if (currentWeight >maxWeight ){
     Info << "current max weight = " << currentWeight << endmsg;
     maxWeight=currentWeight;
   }    
   delete dataEvt;
  }

  hepMCEvtListPtr->rewind();

  Info << "Maximum Weight = " << maxWeight;

  //now do hit and miss method
  EvtSimpleRandomEngine myRandom(3152);
  EvtRandom::setRandomEngine((EvtRandomEngine*)(&myRandom));

  while ( (genEvent = hepMCEvtListPtr->nextEvent()) !=0 ){

    Psi2STo2K2PiGamEvtData* dataEvt=Psi2STo2K2PiGamEvtGenConverter::instance()->convertHepMCtoDataEvt(genEvent);

    _psi2STo2K2PiGamHist->fillHistos("UnselectedHist", dataEvt, 1.); // fill histo for all events before hit & miss
    
    double currentWeight=lhPtr->calcEvtIntensity(dataEvt, theParamVal);
    Info << "current weight = " << currentWeight << endmsg;
    
    double randomNumber=EvtRandom::Flat(maxWeight);
    Info << "random number = " << randomNumber;
    
    if( randomNumber < currentWeight){ //accepted event
      Info << "event accepted!!!" << endmsg;
      genEvent->print(std::cout);
      _selectedGenEvts.push_back(genEvent);
      _psi2STo2K2PiGamHist->fillHistos("SelectedHist", dataEvt, 1.); // fill histo for accepted event          
    }    

    delete dataEvt;
  }
  hepMCEvtListPtr->rewind();

}

Psi2STo2K2PiGamHitAndMiss::~Psi2STo2K2PiGamHitAndMiss()
{
  delete _psi2STo2K2PiGamHist;
}

void Psi2STo2K2PiGamHitAndMiss::dumpToHepMCAscii(const std::string fileName){
 
  HepMC::IO_Ascii ioAscii(fileName.c_str());

  std::vector<HepMC::GenEvent*>::const_iterator it;
  
  for (it=_selectedGenEvts.begin(); it!=_selectedGenEvts.end(); ++it){
    ioAscii.write_event(*it);
  }

}
