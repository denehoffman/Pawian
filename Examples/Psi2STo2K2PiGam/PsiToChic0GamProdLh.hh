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
  PsiToChic0GamProdLh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>);
  PsiToChic0GamProdLh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>);

  /** Destructor */
  virtual ~PsiToChic0GamProdLh();

  virtual AbsPsi2STo2K2PiGamLh* clone_() {
    return new PsiToChic0GamProdLh(_Psi2STo2K2PiGamEvtListPtr);
  }


  // Getters:
  virtual double calcEvtIntensity(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, const param2K2PiGam& theParamVal);
  
  virtual void setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal);
  virtual int setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par);
  virtual unsigned int nFitParams();

  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(param2K2PiGam& theParamVal);
  virtual void dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix);


protected:

  virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);  

  virtual complex<double> chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);  

private:
  unsigned int _nFitParams;
  std::vector<unsigned int> _ampVec;

  complex<double> calcCoherentAmp(Spin Minit, Spin MChi, Spin MGamma, const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);

};

#endif
