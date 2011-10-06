#ifndef _Psi2STo2K2PiGamHitAndMiss_H
#define _Psi2STo2K2PiGamHitAndMiss_H

#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "Examples/Psi2STo2K2PiGam/FitParams2K2PiGam.hh"
#include "HepMC/GenEvent.h"
#include "Utils/PawianCollectionUtils.hh"

#include <map>
#include <boost/shared_ptr.hpp>
#include <vector>




using Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData;

class HepMCEventList;
class AbsPsi2STo2K2PiGamLh;
class Psi2STo2K2PiGamHist;

class Psi2STo2K2PiGamHitAndMiss {

public:

  // create/copy/destroy:

  ///Constructor  
  Psi2STo2K2PiGamHitAndMiss(boost::shared_ptr<AbsPsi2STo2K2PiGamLh> lhPtr, const param2K2PiGam& theParamVal, boost::shared_ptr<HepMCEventList> hepMCEvtListPtr, std::string histName="Psi2SGen.root");



  /** Destructor */
  virtual ~Psi2STo2K2PiGamHitAndMiss();

  // Getters (not working for Psi2STo2K2PiGamHitAndMiss() constructor:
  const std::vector<HepMC::GenEvent*> selectedGenEvts() const {return _selectedGenEvts;}
  void dumpToHepMCAscii(const std::string fileName); 

protected:


private:
  boost::shared_ptr<AbsPsi2STo2K2PiGamLh> _lhPtr;
  std::vector<HepMC::GenEvent*> _selectedGenEvts;
  Psi2STo2K2PiGamHist* _psi2STo2K2PiGamHist;
  

};

#endif 
