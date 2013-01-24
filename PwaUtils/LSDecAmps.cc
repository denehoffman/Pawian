// LSDecAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/LSDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/IsobarLSDecay.hh"
//#include "PwaUtils/XdecAmpRegistry.hh"
#include "Particle/Particle.hh"

#ifdef _OPENMP
#include <omp.h>
#endif

LSDecAmps::LSDecAmps(boost::shared_ptr<IsobarLSDecay> theDec) :
  AbsXdecAmp(theDec)
  ,_JPCLSs(theDec->JPCLSAmps())
  ,_factorMag(1.)
{
  if(_JPCLSs.size()>0) _factorMag=1./sqrt(_JPCLSs.size());
}

LSDecAmps::LSDecAmps(boost::shared_ptr<AbsDecay> theDec) :
  AbsXdecAmp(theDec)
{
}

LSDecAmps::~LSDecAmps()
{
}


complex<double> LSDecAmps::XdecPartAmp(Spin lamX, Spin lamDec, short fixDaughterNr, EvtData* theData, Spin lamFs){

  Spin lam1Min=-_Jdaughter1;
  Spin lam1Max= _Jdaughter1;
  Spin lam2Min=-_Jdaughter2;
  Spin lam2Max=_Jdaughter2;

  if(fixDaughterNr == 1){
     lam1Min = lam1Max = lamDec;
  }
  else if(fixDaughterNr == 2){
     lam2Min = lam2Max = lamDec;
  }
  else{
     Alert << "Invalid fitDaughterNr in XdecPartAmp." << endmsg;
  }

  if( _daughter1IsStable && _Jdaughter1>0){
    lam1Min=lamFs;
    lam1Max=lamFs;
  }
  else if(_daughter2IsStable && _Jdaughter2>0){
    lam2Min=lamFs;
    lam2Max=lamFs;
  }

  complex<double> result(0.,0.);
  std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
  for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){
    if( fabs(lamX) > (*it)->J ) continue;
    double theMag=_currentParamMags[*it];
    double thePhi=_currentParamPhis[*it];
    complex<double> expi(cos(thePhi), sin(thePhi));

    for(Spin lambda1=lam1Min; lambda1<=lam1Max; lambda1++){
      for(Spin lambda2=lam2Min; lambda2<=lam2Max; lambda2++){
        Spin lambda = lambda1-lambda2;
        if( fabs(lambda)>(*it)->J || fabs(lambda)>(*it)->S) continue;

        complex<double> amp = theMag*expi*sqrt(2*(*it)->L+1)
           *Clebsch((*it)->L, 0, (*it)->S, lambda, (*it)->J, lambda)
           *Clebsch(_Jdaughter1, lambda1, _Jdaughter2, -lambda2, (*it)->S, lambda  )
           *conj( theData->WignerDsString[_wignerDKey][(*it)->J][lamX][lambda]);
        result+=amp;
      }
    }
  }

  return result;
}




complex<double> LSDecAmps::XdecAmp(Spin lamX, EvtData* theData, Spin lamFs){
int evtNo=theData->evtNo;
  
  if ( _cacheAmps && !_recalculate){
    complex<double> result(0.,0.);
    result= _cachedAmpMap[evtNo][lamX][lamFs];
    return result;
  }

  Spin lam1Min=-_Jdaughter1;
  Spin lam1Max= _Jdaughter1;
  Spin lam2Min=-_Jdaughter2;
  Spin lam2Max=_Jdaughter2;

  if( _daughter1IsStable && _Jdaughter1>0){
    lam1Min=lamFs;
    lam1Max=lamFs;
  }
  else if(_daughter2IsStable && _Jdaughter2>0){
    lam2Min=lamFs;
    lam2Max=lamFs;
  }

  complex<double> result(0.,0.);
  std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
  for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){
    if( fabs(lamX) > (*it)->J ) continue;
    double theMag=_currentParamMags[*it];
    double thePhi=_currentParamPhis[*it];
    complex<double> expi(cos(thePhi), sin(thePhi));

    for(Spin lambda1=lam1Min; lambda1<=lam1Max; lambda1++){
      for(Spin lambda2=lam2Min; lambda2<=lam2Max; lambda2++){
	Spin lambda = lambda1-lambda2;
	if( fabs(lambda)>(*it)->J || fabs(lambda)>(*it)->S) continue;

	complex<double> amp = theMag*expi*sqrt(2*(*it)->L+1)
	   *Clebsch((*it)->L, 0, (*it)->S, lambda, (*it)->J, lambda)
	   *Clebsch(_Jdaughter1, lambda1, _Jdaughter2, -lambda2, (*it)->S, lambda  )
	   *conj( theData->WignerDsString[_wignerDKey][(*it)->J][lamX][lambda]);

	// std::cout << "amp: " << amp << std::endl;
	amp *=daughterAmp(lambda1, lambda2, theData, lamFs);

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
      _cachedAmpMap[evtNo][lamX][lamFs]=result;
#ifdef _OPENMP
    }
#endif
  }

  return result;
}

void  LSDecAmps::getDefaultParams(fitParams& fitVal, fitParams& fitErr){

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
  for(itLS=_JPCLSs.begin(); itLS!=_JPCLSs.end(); ++itLS){
    currentMagValMap[*itLS]=_factorMag;
    currentPhiValMap[*itLS]=0.;
    currentMagErrMap[*itLS]=0.5*_factorMag;
    currentPhiErrMap[*itLS]=0.3;
  }

  fitVal.Mags[_key]=currentMagValMap;
  fitVal.Phis[_key]=currentPhiValMap;
  fitErr.Mags[_key]=currentMagErrMap;
  fitErr.Phis[_key]=currentPhiErrMap;

  if(_withDyn){
    fitVal.Masses[_massKey]=_decay->motherPart()->mass();
    fitErr.Masses[_massKey]=0.03;
    fitVal.Widths[_massKey]=_decay->motherPart()->width();
    fitErr.Widths[_massKey]=_decay->motherPart()->width();
  }

  if(!_daughter1IsStable) _decAmpDaughter1->getDefaultParams(fitVal, fitErr);
  if(!_daughter2IsStable) _decAmpDaughter2->getDefaultParams(fitVal, fitErr);  
}

void LSDecAmps::print(std::ostream& os) const{
  return; //dummy
}


bool LSDecAmps::checkRecalculation(fitParams& theParamVal){
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
 

void  LSDecAmps::updateFitParams(fitParams& theParamVal){
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
