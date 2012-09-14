#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiToPhiPhiGam/XToPhiPhiDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
// #include "PwaUtils/EvtDataBaseListNew.hh"
#include "Examples/JpsiToPhiPhiGam/JpsiToPhiPhiGamEventList.hh"

XToPhiPhiDecAmps::XToPhiPhiDecAmps(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiToPhiPhiGamStates> theStates, Spin spinX) :
  AbsXdecAmp(name, hypVec, spinX)
  ,_phiPhiKey(name+"ToPhiPhi")
  ,_xBWKey(name+"BreitWigner")
  ,_massIndependent(true)
  ,_theStatesPtr(theStates)
{
  initialize();
}

XToPhiPhiDecAmps::~XToPhiPhiDecAmps()
{
}

complex<double> XToPhiPhiDecAmps::XdecAmp(Spin lamX, EvtDataNew* theData, fitParamsNew& theParamVal){
  complex<double> result(0.,0.);

  result+=XToPhiPhiAmp(lamX, theData, theParamVal);

  complex<double> dynModel(1.,0.);
  if (!_massIndependent){
    double xMass=theParamVal.Masses[_name];
    double xWidth=theParamVal.Widths[_name];
    Vector4<double> p4PhiPhi = theData->FourVecsDec[enumJpsiGamX4V::V4_KsKlKpKm_HeliPsi];
    dynModel=BreitWigner(p4PhiPhi, xMass, xWidth);
  }

  result *=dynModel;
  return result;
}

complex<double> XToPhiPhiDecAmps::XToPhiPhiAmp(Spin lamX, EvtDataNew* theData, fitParamsNew& theParamVal){

  complex<double> result(0.,0.);
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > XToPhiPhiMag=theParamVal.Mags[_phiPhiKey];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > XToPhiPhiPhi=theParamVal.Phis[_phiPhiKey];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;


   for ( itXMag=XToPhiPhiMag.begin(); itXMag!=XToPhiPhiMag.end(); ++itXMag){
     boost::shared_ptr<const JPCLS> XState=itXMag->first;
     double theXMag=itXMag->second;
     double theXPhi=XToPhiPhiPhi[XState];
     complex<double> expiphiX(cos(theXPhi), sin(theXPhi));

    for(Spin lambdaPhi1=-1; lambdaPhi1<=1; lambdaPhi1++){
      for(Spin lambdaPhi2=-1; lambdaPhi2<=1; lambdaPhi2++){
        Spin lambda = lambdaPhi1-lambdaPhi2;
        if( fabs(lambda)>XState->J || fabs(lambda)>XState->S) continue;
        
        complex<double> amp = theXMag*expiphiX*sqrt(2*XState->L+1)
          *Clebsch(XState->L, 0, XState->S, lambda, XState->J, lambda)
          *Clebsch(1, lambdaPhi1, 1, -lambdaPhi2, XState->S, lambda  )
          *conj( theData->WignerDsDec[enumJpsiGamXDfunc::Df_XToPhiPhi][XState->J][_J_X][lambda]  );

        
        amp = amp * phiphiTo4KAmp( theData, lambdaPhi1, lambdaPhi2 );

        result +=amp;
      }
    }
   }
   
   result*=conj(theData->WignerDsDec[enumJpsiGamXDfunc::Df_XToPhiPhi][_J_X][lamX][0]);
   return result;
}

complex<double> XToPhiPhiDecAmps::phiphiTo4KAmp( EvtDataNew* theData, Spin lambdaPhi1, Spin lambdaPhi2 ){
  complex<double> result(0.,0.);
  
  result = 3. * conj(theData->WignerDsDec[enumJpsiGamXDfunc::Df_PhiToKsKl][1][lambdaPhi1][0])
    * 3.* conj(theData->WignerDsDec[enumJpsiGamXDfunc::Df_PhiToKpKm][1][lambdaPhi2][0]);
  
  return result;
}


void  XToPhiPhiDecAmps::getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr){

  std::vector< boost::shared_ptr<const JPCLS> > PhiPhiStates;
  if(_J_X==0) PhiPhiStates=_theStatesPtr->EtaToPhiPhiStates();
  else if(_J_X==1) PhiPhiStates=_theStatesPtr->F1ToPhiPhiStates();
  else if(_J_X==2) PhiPhiStates=_theStatesPtr->Eta2ToPhiPhiStates();  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
  for(itLS=PhiPhiStates.begin(); itLS!=PhiPhiStates.end(); ++itLS){
    currentMagValMap[*itLS]=0.2;
    currentPhiValMap[*itLS]=0.;
    currentMagErrMap[*itLS]=0.8;
    currentPhiErrMap[*itLS]=0.3;
  }
  fitVal.Mags[_phiPhiKey]=currentMagValMap;
  fitVal.Phis[_phiPhiKey]=currentPhiValMap;
  fitErr.Mags[_phiPhiKey]=currentMagErrMap;
  fitErr.Phis[_phiPhiKey]=currentPhiErrMap;
  
  
  if (!_massIndependent){
    size_t pos=_name.find("_");
    std::string massMeVString=_name.substr(pos+1); 
    stringstream massMeVStrStream(massMeVString);
    int MassMeV;
    massMeVStrStream >> MassMeV;
    double MassGeV= ( (double) MassMeV)/1000.;
    
    fitVal.Masses[_name]=MassGeV;
    fitErr.Masses[_name]=0.01;
    fitVal.Widths[_name]=0.2;
    fitErr.Widths[_name]=0.02;
  }
}

void XToPhiPhiDecAmps::print(std::ostream& os) const{
  return; //dummy
}

void XToPhiPhiDecAmps::initialize(){
  std::vector<std::string>::const_iterator it;

  for (it=_hypVec.begin(); it!=_hypVec.end(); ++it){

    if (it->compare(0, _xBWKey.size(), _xBWKey) ==0){
      _massIndependent=false;
    }
    
  }

}
