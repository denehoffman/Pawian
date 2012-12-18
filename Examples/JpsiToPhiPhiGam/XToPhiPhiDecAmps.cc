#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiToPhiPhiGam/XToPhiPhiDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
// #include "PwaUtils/EvtDataBaseListNew.hh"
#include "Examples/JpsiToPhiPhiGam/JpsiToPhiPhiGamEventList.hh"

#ifdef _OPENMP
#include <omp.h>
#endif

XToPhiPhiDecAmps::XToPhiPhiDecAmps(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiToPhiPhiGamStates> theStates, Spin spinX, int parity) :
  AbsXdecAmp(name, hypVec, spinX, parity)
  ,_phiPhiKey(name+"ToPhiPhi")
  ,_xBWKey(name+"BreitWigner")
  ,_xFlatteKey(name+"Flatte")
  ,_massIndependent(true)
  ,_bwMassFit(false)
  ,_flatteMassFit(false)
  ,_phiMass( 1.019455)
  ,_omegaMass( 0.78265)
  ,_phiPhiPair(make_pair(_phiMass, _phiMass))
  ,_omegaPhiPair(make_pair(_omegaMass, _phiMass))
  ,_theStatesPtr(theStates)
{
  initialize();
}

XToPhiPhiDecAmps::~XToPhiPhiDecAmps()
{
}

complex<double> XToPhiPhiDecAmps::XdecAmp(Spin lamX, EvtDataNew* theData){
  int evtNo=theData->evtNo;
  
  if ( _cacheAmps && !_recalculate) return _cachedAmpMap[evtNo][lamX];


  complex<double> result(0.,0.);

    result+=XToPhiPhiAmp(lamX, theData);
    
    complex<double> dynModel(1.,0.);
    if (!_massIndependent){
      Vector4<double> p4PhiPhi = theData->FourVecsDec[enumJpsiGamX4V::V4_KsKlKpKm_HeliPsi];
      if(_bwMassFit){
	dynModel=BreitWigner(p4PhiPhi, _currentXMass, _currentXWidth);
      }
      else if(_flatteMassFit){
	dynModel=Flatte( p4PhiPhi, _phiPhiPair, _omegaPhiPair, _currentXMass, _currentgFactorPhiPhi, _currentgFactorOmegaPhi  );
      }
    }
    
    result *=dynModel;

    if ( _cacheAmps){    
#ifdef _OPENMP
#pragma omp critical
      {
#endif
	_cachedAmpMap[evtNo][lamX]=result;
#ifdef _OPENMP
      }
#endif
    }
  
  return result;
}

complex<double> XToPhiPhiDecAmps::XToPhiPhiAmp(Spin lamX, EvtDataNew* theData){

  complex<double> result(0.,0.);
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;


   for ( itXMag=_currentParamMags.begin(); itXMag!=_currentParamMags.end(); ++itXMag){
     boost::shared_ptr<const JPCLS> XState=itXMag->first;
     double theXMag=itXMag->second;
     double theXPhi=_currentParamPhis[XState];
     complex<double> expiphiX(cos(theXPhi), sin(theXPhi));

    for(Spin lambdaPhi1=-1; lambdaPhi1<=1; lambdaPhi1++){
      for(Spin lambdaPhi2=-1; lambdaPhi2<=1; lambdaPhi2++){
        Spin lambda = lambdaPhi1-lambdaPhi2;
        if( fabs(lambda)>XState->J || fabs(lambda)>XState->S) continue;
        
        complex<double> amp = theXMag*expiphiX*sqrt(2*XState->L+1)
          *Clebsch(XState->L, 0, XState->S, lambda, XState->J, lambda)
          *Clebsch(1, lambdaPhi1, 1, -lambdaPhi2, XState->S, lambda  )
          *conj( theData->WignerDsDec[enumJpsiGamXDfunc::Df_XToPhiPhi1][XState->J][lamX][lambda]);
        
        amp = amp * phiphiTo4KAmp( theData, lambdaPhi1, lambdaPhi2 );

        result +=amp;
      }
    }
   }
 
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
  if(_J_X==0 && _parity==-1) PhiPhiStates=_theStatesPtr->EtaToPhiPhiStates();
  else if(_J_X==0 && _parity==1) PhiPhiStates=_theStatesPtr->F0ToPhiPhiStates();
  else if(_J_X==1 && _parity==1) PhiPhiStates=_theStatesPtr->F1ToPhiPhiStates();
  else if(_J_X==2 && _parity==-1) PhiPhiStates=_theStatesPtr->Eta2ToPhiPhiStates();
  else if(_J_X==2 && _parity==1) PhiPhiStates=_theStatesPtr->F2ToPhiPhiStates();  
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

    if(_bwMassFit){
      fitVal.Widths[_name]=0.2;
      fitErr.Widths[_name]=0.02;
    }
    else if( _flatteMassFit){
      std::string phiPhiGFactorKey="gFactorPhiPhi_"+_name;
      std::string omegaPhiGFactorKey="gFactorOmegaPhi_"+_name;
      fitVal.gFactors[phiPhiGFactorKey]=0.3;
      fitVal.gFactors[omegaPhiGFactorKey]=0.3;
      fitErr.gFactors[phiPhiGFactorKey]=0.2;
      fitErr.gFactors[omegaPhiGFactorKey]=0.2;
    }
  }
}

void XToPhiPhiDecAmps::print(std::ostream& os) const{
  return; //dummy
}

void XToPhiPhiDecAmps::initialize(){
  std::vector<std::string>::const_iterator it;

  for (it=_hypVec.begin(); it!=_hypVec.end(); ++it){
    if (it->compare(0, _xBWKey.size(), _xBWKey) ==0){
      _bwMassFit=true;
      _massIndependent=false;
    }
    else if (it->compare(0, _xFlatteKey.size(), _xFlatteKey) ==0){
      _flatteMassFit=true;
      _massIndependent=false;
    }    
  }


}

bool XToPhiPhiDecAmps::checkRecalculation(fitParamsNew& theParamVal){
  _recalculate=false;
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > XToPhiPhiMag=theParamVal.Mags[_phiPhiKey];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > XToPhiPhiPhi=theParamVal.Phis[_phiPhiKey];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;
   for ( itXMag=XToPhiPhiMag.begin(); itXMag!=XToPhiPhiMag.end(); ++itXMag){
     boost::shared_ptr<const JPCLS> XState=itXMag->first;
     double theXMag=itXMag->second;
     double theXPhi=XToPhiPhiPhi[XState];
     if ( fabs(theXMag - _currentParamMags[XState])  > 1.e-10 ){
       _recalculate=true;
     }
     if ( fabs(theXPhi - _currentParamPhis[XState])  > 1.e-10 ){
       _recalculate=true;
     }
   }

   if (!_massIndependent){
     double xMass=theParamVal.Masses[_name];
     if ( fabs(xMass-_currentXMass) > 1.e-10){
       _recalculate=true;
     } 
    
     if(_bwMassFit){
       double xWidth=theParamVal.Widths[_name];
       if ( fabs(xWidth - _currentXWidth) > 1.e-10){
	 _recalculate=true;
       }      
     }
    else if(_flatteMassFit){
      std::string phiPhiGFactorKey="gFactorPhiPhi_"+_name;
      std::string omegaPhiFactorKey="gFactorOmegaPhi_"+_name;
      double gFactorPhiPhi=theParamVal.gFactors[phiPhiGFactorKey];
      double gFactorOmegaPhi=theParamVal.gFactors[omegaPhiFactorKey];
      if ( fabs(gFactorPhiPhi - _currentgFactorPhiPhi)  > 1.e-10 ){
	_recalculate=true;
      }
      if ( fabs(gFactorOmegaPhi - _currentgFactorOmegaPhi)  > 1.e-10 ){
	_recalculate=true;
      }
    }
  }
   if (_recalculate) Info << "Recalculate amplitude:\t" << _name << endmsg;  

   return _recalculate;  
}
 

void  XToPhiPhiDecAmps::updateFitParams(fitParamsNew& theParamVal){
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > XToPhiPhiMag=theParamVal.Mags[_phiPhiKey];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > XToPhiPhiPhi=theParamVal.Phis[_phiPhiKey];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;
   for ( itXMag=XToPhiPhiMag.begin(); itXMag!=XToPhiPhiMag.end(); ++itXMag){
     boost::shared_ptr<const JPCLS> XState=itXMag->first;
     double theXMag=itXMag->second;
     double theXPhi=XToPhiPhiPhi[XState];
     _currentParamMags[XState]=theXMag;
     _currentParamPhis[XState]=theXPhi;
   }

   if (!_massIndependent){
     double xMass=theParamVal.Masses[_name];
     _currentXMass= xMass; 
    
     if(_bwMassFit){
       double xWidth=theParamVal.Widths[_name];
       _currentXWidth=xWidth;      
     }
    else if(_flatteMassFit){
      std::string phiPhiGFactorKey="gFactorPhiPhi_"+_name;
      std::string omegaPhiFactorKey="gFactorOmegaPhi_"+_name;
      double gFactorPhiPhi=theParamVal.gFactors[phiPhiGFactorKey];
      double gFactorOmegaPhi=theParamVal.gFactors[omegaPhiFactorKey];
      _currentgFactorPhiPhi=gFactorPhiPhi;
      _currentgFactorOmegaPhi=gFactorOmegaPhi;
    }
  }
}
