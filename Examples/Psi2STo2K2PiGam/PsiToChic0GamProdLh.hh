#ifndef _PsiToChic0GamProdLh_H
#define _PsiToChic0GamProdLh_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "Examples/Psi2STo2K2PiGam/AbsPsi2STo2K2PiGamLh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"



class PsiToChic0GamProdLh : public AbsPsi2STo2K2PiGamLh{

public:

  // create/copy/destroy:

  ///Constructor 
  PsiToChic0GamProdLh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>, boost::shared_ptr<const Psi2STo2K2PiGamStates>);
  PsiToChic0GamProdLh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>);

  /** Destructor */
  virtual ~PsiToChic0GamProdLh();

  virtual AbsPsi2STo2K2PiGamLh* clone_() const {
    return new PsiToChic0GamProdLh(_Psi2STo2K2PiGamEvtListPtr, _Psi2STo2K2PiGamStatesPtr);
  }


  // Getters:
  virtual double calcEvtIntensity(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, const Psi2STo2K2PiGamData::fitParamVal& theParamVal);
  
  virtual void setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal);
  virtual int setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const;

  virtual void print(std::ostream& os) const;
  //  virtual void printCurrentFitResult(Psi2STo2K2PiGamData::fitParamVal& theParamVal) const;

protected:

  virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);  virtual complex<double> chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);  

private:

  complex<double> calcCoherentAmp(Spin Minit, Spin MChi, Spin MGamma, const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);

};

#endif
