#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamKsKlKK/AbsJpsiGamKsKlKKLh.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKEventList.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

AbsJpsiGamKsKlKKLh::AbsJpsiGamKsKlKKLh(boost::shared_ptr<const JpsiGamKsKlKKEventList> theEvtList) :
  _JpsiGamKsKlKKEvtListPtr(theEvtList),
  _fitparamsGamKsKlKK()
{
  _evtDataVec=_JpsiGamKsKlKKEvtListPtr->getDataVecs();
  _evtMCVec=_JpsiGamKsKlKKEvtListPtr->getMcVecs();
}

AbsJpsiGamKsKlKKLh::AbsJpsiGamKsKlKKLh(boost::shared_ptr<AbsJpsiGamKsKlKKLh> theAbsJpsiGamKsKlKKLhPtr):
  _JpsiGamKsKlKKEvtListPtr(theAbsJpsiGamKsKlKKLhPtr->getEventList()),
  _fitparamsGamKsKlKK()
{
  _evtDataVec=_JpsiGamKsKlKKEvtListPtr->getDataVecs();
  _evtMCVec=_JpsiGamKsKlKKEvtListPtr->getMcVecs();
}

AbsJpsiGamKsKlKKLh::~AbsJpsiGamKsKlKKLh()
{
}

double AbsJpsiGamKsKlKKLh::calcLogLh(const paramGamKsKlKK& theParamVal){
 
  double logLH=0.;
  double logLH_data=0.;

  std::vector<JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData*>::iterator iterd;
  for (iterd=_evtDataVec.begin(); iterd!=_evtDataVec.end(); ++iterd){
    double intensity=calcEvtIntensity((*iterd), theParamVal);
    if (intensity>0.) logLH_data+=log10(intensity);
  } 

  double LH_mc=0.;
  
  std::vector<JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData*>::iterator iterm;
  for (iterm=_evtMCVec.begin(); iterm!=_evtMCVec.end(); ++iterm){
           double intensity=calcEvtIntensity((*iterm), theParamVal);
           LH_mc+=intensity;
         }

  double logLH_mc_Norm=0.;
  if (LH_mc>0.) logLH_mc_Norm=log10(LH_mc/_evtMCVec.size());

  logLH=_evtDataVec.size()/2.*(LH_mc/_evtMCVec.size()-1)*(LH_mc/_evtMCVec.size()-1)
    -logLH_data
    +_evtDataVec.size()*logLH_mc_Norm;

  Info << "current LH = " << logLH ;  // << endmsg;

 return logLH;

}

double AbsJpsiGamKsKlKKLh::calcEvtIntensity(JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData, const paramGamKsKlKK& theParamVal){

  double phaseSpaceVal=0.0 ; //theParamVal.phaseSpace;
  
  Spin Psi2SM=1;
  Spin GamM=1;
  complex<double> AmpPsi2SMpGp=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData);

  Psi2SM=1;
  GamM=-1;
  complex<double> AmpPsi2SMpGm=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData);

  Psi2SM=-1;
  GamM=1; 
  complex<double> AmpPsi2SMmGp=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData);

  Psi2SM=-1;
  GamM=-1; 
  complex<double> AmpPsi2SMmGm=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData);

//   DebugMsg << "AmpPsi2SMp " << AmpPsi2SMp ;  // << endmsg;
  
  double result=norm(AmpPsi2SMpGp)+norm(AmpPsi2SMpGm)+norm(AmpPsi2SMmGp)+norm(AmpPsi2SMmGm)+phaseSpaceVal;

  return result;  
}




void AbsJpsiGamKsKlKKLh::print(std::ostream& os) const{
  os << "AbsJpsiGamKsKlKKLh::print\n";
}



