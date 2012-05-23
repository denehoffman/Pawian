#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamEtaPiPiNew/F1DecAmp.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiData.hh"

F1DecAmp::F1DecAmp(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates) :
  XDecAmpBase(name, hypVec, theStates, Spin(1))
{
}

F1DecAmp::~F1DecAmp()
{
}

complex<double> F1DecAmp::XdecAmp(Spin lamX, EvtData* theData, fitParamsNew& theParamVal){
  complex<double> result(0.,0.);

  if(_a980piHyp){
    result+=XToAPiFlatteAmp(lamX, theData, theParamVal);
  }
  if(_sigmaEtaHyp){
    double sigmaMass=theParamVal.Masses["sigma"];
    double sigmaWidth=theParamVal.Widths["sigma"];
    result+=XToEtaFAmp(lamX, 1, theData, theParamVal.Mags[_sigmaEtaKey], theParamVal.Phis[_sigmaEtaKey],sigmaMass, sigmaWidth);
  }

  return result;
}


void  F1DecAmp::getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr){

  if(_a980piHyp){

    std::vector< boost::shared_ptr<const JPCLS> > a0PiStates=_theStatesPtr->F1Toa0PiStates();

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
    
    for(itLS=a0PiStates.begin(); itLS!=a0PiStates.end(); ++itLS){
      currentMagValMap[*itLS]=0.2;
      currentPhiValMap[*itLS]=0.;
      currentMagErrMap[*itLS]=0.8;
      currentPhiErrMap[*itLS]=0.3;
    }
    fitVal.Mags[_a980piKey]=currentMagValMap;
    fitVal.Phis[_a980piKey]=currentPhiValMap;
    fitErr.Mags[_a980piKey]=currentMagErrMap;
    fitErr.Phis[_a980piKey]=currentPhiErrMap;

    fitVal.Masses["a0_980"]=0.98;
    fitErr.Masses["a0_980"]=0.03;
    fitVal.gFactors["a0_980gPiEta"]=1.;
    fitErr.gFactors["a0_980gPiEta"]=0.2;
    fitVal.gFactors["a0_980gKK"]=1.;
    fitErr.gFactors["a0_980gKK"]=0.2;
  }

  if(_sigmaEtaHyp){
    std::vector< boost::shared_ptr<const JPCLS> > f1EtaStates=_theStatesPtr->F1Tof0EtaStates();
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
    
    for(itLS=f1EtaStates.begin(); itLS!=f1EtaStates.end(); ++itLS){
      currentMagValMap[*itLS]=0.2;
      currentPhiValMap[*itLS]=0.;
      currentMagErrMap[*itLS]=0.8;
      currentPhiErrMap[*itLS]=0.3;
    }
    fitVal.Mags[_sigmaEtaKey]=currentMagValMap;
    fitVal.Phis[_sigmaEtaKey]=currentPhiValMap;
    fitErr.Mags[_sigmaEtaKey]=currentMagErrMap;
    fitErr.Phis[_sigmaEtaKey]=currentPhiErrMap;

    fitVal.Masses["sigma"]=0.6;
    fitErr.Masses["sigma"]=0.1;
    fitVal.Widths["sigma"]=0.5;
    fitErr.Widths["sigma"]=0.2;

  }
}

void F1DecAmp::print(std::ostream& os) const{
  return; //dummy
}


