#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarp/pbarpDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/DataUtils.hh"
#include "Examples/pbarp/IsobarDecay.hh"
#include "Particle/Particle.hh"

#ifdef _OPENMP
#include <omp.h>
#endif

pbarpDecAmps::pbarpDecAmps(boost::shared_ptr<IsobarDecay> theDec) :
  AbsXdecAmp(theDec->name())
  ,_decay(theDec)
  ,_JPCPtr(theDec->motherJPC())
  ,_JPCLSs(theDec->JPCLSAmps())
  ,_key("_"+theDec->fitParSuffix())
  ,_massKey("")
  ,_wignerDKey(theDec->wignerDKey())
  ,_daughter1IsStable(theDec->isDaughter1Stable())
  ,_daughter2IsStable(theDec->isDaughter2Stable())
  ,_withDyn(theDec->withDynamics())
{
  initialize();
}

pbarpDecAmps::~pbarpDecAmps()
{
}

complex<double> pbarpDecAmps::XdecAmp(Spin lamX, EvtDataNew* theData){
int evtNo=theData->evtNo;
  
 if ( _cacheAmps && !_recalculate){
   complex<double> result(0.,0.);    
#ifdef _OPENMP
#pragma omp critical
   {
#endif
     result= _cachedAmpMap[evtNo][lamX];
#ifdef _OPENMP
     
   }
#endif
   return result;
 }

  complex<double> result(0.,0.);
   std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
   for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){
     if( fabs(lamX) > (*it)->J ) continue;
     double theMag=_currentParamMags[*it];
     double thePhi=_currentParamPhis[*it];
     complex<double> expi(cos(thePhi), sin(thePhi));

    for(Spin lambda1=-_Jdaughter1; lambda1<=_Jdaughter1; lambda1++){
      for(Spin lambda2=-_Jdaughter2; lambda2<=_Jdaughter2; lambda2++){
	Spin lambda = lambda1-lambda2;
	if( fabs(lambda)>(*it)->J || fabs(lambda)>(*it)->S) continue;

	complex<double> amp = theMag*expi*sqrt(2*(*it)->L+1)
          *Clebsch((*it)->L, 0, (*it)->S, lambda, (*it)->J, lambda)
          *Clebsch(_Jdaughter1, lambda1, _Jdaughter2, -lambda2, (*it)->S, lambda  )
          *conj( theData->WignerDsString[_wignerDKey][(*it)->J][lamX][lambda]);

	// std::cout << "amp: " << amp << std::endl;
	amp *=daughterAmp(lambda1, lambda2, theData);

	result+=amp;
      }
    }
   }

   if(_withDyn){
     Vector4<double> mass4Vec(0.,0.,0.,0.);
     std::vector<Particle*> fsParticleVec=_decay->finalStateParticles();

     std::vector<Particle*>::iterator itPartVec;
     for (itPartVec=fsParticleVec.begin(); itPartVec!=fsParticleVec.end(); ++itPartVec){
       mass4Vec+=theData->FourVecsString[(*itPartVec)->name()];
     }
     result*=BreitWigner(mass4Vec, _currentXMass, _currentXWidth);     
   }
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

complex<double> pbarpDecAmps::daughterAmp(Spin lam1, Spin lam2, EvtDataNew* theData){
  complex<double> result(1.,0.);
  if(!_daughter1IsStable) result *= _decAmpDaughter1->XdecAmp(lam1, theData);
  if(!_daughter2IsStable) result *= _decAmpDaughter2->XdecAmp(lam2, theData);
  return result;
}


void  pbarpDecAmps::getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr){

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
  for(itLS=_JPCLSs.begin(); itLS!=_JPCLSs.end(); ++itLS){
    currentMagValMap[*itLS]=1.;
    currentPhiValMap[*itLS]=0.;
    currentMagErrMap[*itLS]=0.5;
    currentPhiErrMap[*itLS]=0.3;
  }

  fitVal.Mags[_key]=currentMagValMap;
  fitVal.Phis[_key]=currentPhiValMap;
  fitErr.Mags[_key]=currentMagErrMap;
  fitErr.Phis[_key]=currentPhiErrMap;

  if(_withDyn){
    fitVal.Masses[_massKey]=_decay->motherPart()->mass();
    fitErr.Masses[_massKey]=_decay->motherPart()->width();
    fitVal.Widths[_massKey]=_decay->motherPart()->width();
    fitErr.Widths[_massKey]=_decay->motherPart()->width();
  }

  if(!_daughter1IsStable) _decAmpDaughter1->getDefaultParams(fitVal, fitErr);
  if(!_daughter2IsStable) _decAmpDaughter2->getDefaultParams(fitVal, fitErr);  
}

void pbarpDecAmps::print(std::ostream& os) const{
  return; //dummy
}

void pbarpDecAmps::initialize(){

  if(_withDyn){
    if(!_decay->hasMother()){
      Alert << "no mother resonance; can not add dynamis" << endmsg;
      exit(1);
    }
    _massKey=_decay->motherPart()->name();
  }
  
  if(!_daughter1IsStable){
    boost::shared_ptr<IsobarDecay> decDaughter1=_decay->decDaughter1();
    _decAmpDaughter1=boost::shared_ptr<pbarpDecAmps>(new pbarpDecAmps(decDaughter1));
  }
  
  if(!_daughter2IsStable){
    boost::shared_ptr<IsobarDecay> decDaughter2=_decay->decDaughter1();
    _decAmpDaughter2=boost::shared_ptr<pbarpDecAmps>(new pbarpDecAmps(decDaughter2));
  }
  
  _Jdaughter1=(Spin) _decay->daughter1Part()->J();
  _Jdaughter2=(Spin) _decay->daughter2Part()->J();
}

bool pbarpDecAmps::checkRecalculation(fitParamsNew& theParamVal){
  _recalculate=false;

   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > magMap=theParamVal.Mags[_key];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > phiMap=theParamVal.Phis[_key];

   std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
   for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){
     double theMag=magMap[*it];
     double thePhi=phiMap[*it];

     if ( fabs(theMag - _currentParamMags[*it])  > 1.e-10 ){
       _recalculate=true;
     }
     if ( fabs(thePhi - _currentParamPhis[*it])  > 1.e-10 ){
       _recalculate=true;
     }
   }

   if(_withDyn){
     double mass=theParamVal.Masses[_massKey];
     if ( fabs(mass-_currentXMass) > 1.e-10){
       _recalculate=true;
     }
     double width=theParamVal.Widths[_massKey];
     if ( fabs(width-_currentXWidth) > 1.e-10){
       _recalculate=true;
     }
   }

   if(!_daughter1IsStable) {
     if(_decAmpDaughter1->checkRecalculation(theParamVal)) _recalculate=true;
   }
   if(!_daughter2IsStable){
     if(_decAmpDaughter2->checkRecalculation(theParamVal)) _recalculate=true;
   }
   return _recalculate;
}
 

void  pbarpDecAmps::updateFitParams(fitParamsNew& theParamVal){
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > magMap=theParamVal.Mags[_key];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > phiMap=theParamVal.Phis[_key];

   std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
   for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){
     double theMag=magMap[*it];
     double thePhi=phiMap[*it];
     _currentParamMags[*it]=theMag;
     _currentParamPhis[*it]=thePhi;
   }

   if (_withDyn){
     double xMass=theParamVal.Masses[_massKey];
     _currentXMass= xMass;
     double xWidth=theParamVal.Widths[_massKey];
     _currentXWidth=xWidth;
   }

  if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(theParamVal);
  if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(theParamVal);

}
