//************************************************************************//
//                                                                        //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)                      //
//               Julian Pychy (julian@ep1.rub.de)                         //
//               - Ruhr-Universität Bochum                                //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

// gammapBaseLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <fstream>
#include <getopt.h>
#include <math.h>
#include <set>
#include <string>

#include "ErrLogger/ErrLogger.hh"
#include "ConfigParser/ParserBase.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "Particle/Particle.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/FsParticleProjections.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "gammapUtils/GammapChannelEnv.hh"
#include "gammapUtils/ReflectivityDecAmps.hh"
#include "gammapUtils/gammapBaseLh.hh"
#include "gammapUtils/gammapReaction.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

gammapBaseLh::gammapBaseLh(ChannelID channelID)
    : AbsLh(channelID), _useReflectivity(false) {
  initialize();
}

gammapBaseLh::~gammapBaseLh() { ; }

complex<double> gammapBaseLh::calcSpinDensity(Spin M1, Spin M2,
                                              std::string &nameDec,
                                              EvtData *theData) {

  complex<double> result(0., 0.);

  return result;
}

complex<double>
gammapBaseLh::calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec,
                              EvtData *theData,
                              std::map<std::shared_ptr<const JPCLS>,
                                       std::vector<std::shared_ptr<AbsXdecAmp>>,
                                       pawian::Collection::SharedPtrLess>
                                  pbarpAmps) {
  complex<double> resultAmp(0., 0.);

  return resultAmp;
}

double
gammapBaseLh::calcEvtIntensity(EvtData *theData,
                               std::shared_ptr<AbsPawianParameters> fitPar) {
  double result = 0.;

  std::vector<std::shared_ptr<AbsXdecAmp>>::iterator itDecAll;
  for (itDecAll = _decAmps.begin(); itDecAll != _decAmps.end(); ++itDecAll) {
    (*itDecAll)->calcDynamics(theData);
  }

  const std::vector<std::vector<Spin>> &spinProjections =
      _fsParticleProjections->spinProjections();
  const std::array<Spin, 2> photonHelicities = {Spin(1), Spin(-1)};
  const std::array<Spin, 2> targetSpinProjections = {Spin(0.5), Spin(-0.5)};

  for (unsigned int projId = 0; projId < spinProjections.size(); ++projId) {
    const std::vector<Spin> &currentSpinProjection = spinProjections.at(projId);
    for (itDecAll = _decAmps.begin(); itDecAll != _decAmps.end(); ++itDecAll) {
      (*itDecAll)->setSpinProjections(currentSpinProjection, projId);
    }

    for (unsigned int targetId = 0; targetId < targetSpinProjections.size();
         ++targetId) {
      std::array<std::complex<double>, 2> helicityAmps = {
          std::complex<double>(0., 0.), std::complex<double>(0., 0.)};

      for (unsigned int photonId = 0; photonId < photonHelicities.size();
           ++photonId) {
        if (_useReflectivity) {
          // The target moves opposite to the beam, so its physical helicity
          // is the negative of its spin projection on the beam axis.
          const Spin initialProtonHelicity =
              -targetSpinProjections.at(targetId);
          for (itDecAll = _decAmps.begin(); itDecAll != _decAmps.end();
               ++itDecAll) {
            std::shared_ptr<ReflectivityDecAmps> reflectivityAmp =
                std::dynamic_pointer_cast<ReflectivityDecAmps>(*itDecAll);
            if (!reflectivityAmp) {
              Alert << "Non-reflectivity production amplitude in a "
                       "Reflectivity channel"
                    << endmsg;
              exit(1);
            }
            helicityAmps.at(photonId) += reflectivityAmp->amplitude(
                photonHelicities.at(photonId), initialProtonHelicity,
                theData);
          }
        } else {
          const Spin lamGammap = photonHelicities.at(photonId) +
                                 targetSpinProjections.at(targetId);
          for (itDecAll = _decAmps.begin(); itDecAll != _decAmps.end();
               ++itDecAll) {
            helicityAmps.at(photonId) +=
                (*itDecAll)->XdecAmp(lamGammap, theData);
          }
        }
      }

      result += _beamPolarization.intensity(helicityAmps);
    }
  }

  if (_usePhasespace)
    result += fitPar->Value(_phasespaceKey);

  result *= fitPar->Value(_channelScaleParam);
  return result;
}

// double gammapBaseLh::calcEvtIntensity(EvtData* theData, fitParCol&
// theParamVal){

//   double result=0.;
//   std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
//   std::map< Spin, std::map< Spin, std::map< Spin, complex<double> > > >
//   mapPinGamiPout; for (Spin mzPin=-1./2; mzPin<=1./2; mzPin++){ 	for
//   (Spin mzGam=-1.; mzGam<=1.; mzGam+=2){ 	  for (Spin lamPout=-1./2;
//   lamPout<=1./2; lamPout++){
//   mapPinGamiPout[mzPin][mzGam][lamPout]=complex<double> (0.,0.);
//   	  }
//   	}
//   }

//   std::map< std::shared_ptr<const jpcRes>, std::vector< std::shared_ptr<const
//   JPCLS> >, pawian::Collection::SharedPtrLess>::iterator itjpc; std::vector<
//   std::shared_ptr<const JPCLS> >::iterator itJPClj;

//   for(itDec=_decAmps.begin(); itDec !=_decAmps.end(); ++itDec){
//     std::shared_ptr<const jpcRes> currentJPC=(*itDec)->jpcPtr();
//     double iso12Val=_currentParamJPCIsos12.at(currentJPC);
//     double iso32Val=_currentParamJPCIsos32.at(currentJPC);
//     double isoFactor=iso12Val;
//     if( int(2.*(*itDec)->absDec()->motherIGJPC()->I)==3) isoFactor=iso32Val;

//     double theMag=_currentParamMags.at(currentJPC);
//     double thePhi=0.;
//     if( 2.*currentJPC->J > 1.) thePhi=_currentParamPhis.at(currentJPC);
//     else if(currentJPC->P == 1){ // 1/2+ M1
//       thePhi=M_PI/2.;
//     }

//     complex<double> current_amp_prod(0.,0.);
//     std::vector< std::shared_ptr<const JPCLS> > currentJPCljVec =
//     _jpcToJPCljMap.at(currentJPC); for(itJPClj=currentJPCljVec.begin();
//     itJPClj!=currentJPCljVec.end(); ++itJPClj){
//       int lmp=(*itJPClj)->L;
//       Spin jmp=(*itJPClj)->S;
//       Spin J=(*itJPClj)->J;

//       double elMagFactor=1.;
//       if(currentJPC->P == pow(-1.,(int)(J + 1./2))){ // 1/2-,3/2+,5/2-,...
// 	if(jmp==J+Spin(1./2)){ //el. multipole
// 	  elMagFactor=cos(thePhi);
// 	}
// 	else elMagFactor=sin(thePhi); //mag multipole
//       }
//       else{ // 1/2+,3/2-,5/2+,...
// 	if(jmp==J+Spin(1./2)){ //mag. multipole
// 	  elMagFactor=sin(thePhi);
// 	}
// 	else elMagFactor=cos(thePhi);  //el. multipole
//       }

//       for (Spin mzPin=-1./2; mzPin<=1./2; mzPin++){
// 	for (Spin mzGam=-1.; mzGam<=1.; mzGam+=2.){
//   	  Spin mz = mzGam + mzPin; // spin-j projection
// 	  if( J < fabs(mz) || jmp < fabs(mzGam) ) continue;
//   	  complex<double>
//   current_amp_prod=theMag*elMagFactor*isoFactor*sqrt(2*lmp+1)*Clebsch(1,
//   mzGam,lmp,0,jmp,mzGam)*Clebsch(jmp, mzGam,1/2., mzPin,J,mz);
// 	  for (Spin lamPout=-1./2; lamPout<=1./2; lamPout++){
//   	    mapPinGamiPout.at(mzPin).at(mzGam).at(lamPout)+=current_amp_prod*(*itDec)->XdecAmp(mz,
//   theData, lamPout);
//   	  }
//   	}
//       }
//     }
//   }

//   std::map< Spin, std::map< Spin, std::map< Spin, complex<double> > >
//   >::iterator itAmpMap;
//     for(itAmpMap=mapPinGamiPout.begin(); itAmpMap!=mapPinGamiPout.end();
//     ++itAmpMap){
//       std::map< Spin, std::map< Spin, complex<double> > >::iterator
//       itAmpMap2; for(itAmpMap2=itAmpMap->second.begin();
//       itAmpMap2!=itAmpMap->second.end(); ++itAmpMap2){
//   	std::map< Spin, complex<double> >::iterator itAmpMap3;
//   	for(itAmpMap3=itAmpMap2->second.begin();
//   itAmpMap3!=itAmpMap2->second.end(); ++itAmpMap3){
//   	  result+=norm(itAmpMap3->second);
//   	}
//       }
//   }

//   if(_usePhasespace) result+=theParamVal.otherParams[_phasespaceKey];

//   result *= theParamVal.otherParams.at(_channelScaleParam);

//   return result;

// }

// void gammapBaseLh::getDefaultParams(fitParCol& fitVal, fitParCol& fitErr){
//   AbsLh::getDefaultParams(fitVal, fitErr);

//   std::map< std::shared_ptr<const jpcRes>, double,
//   pawian::Collection::SharedPtrLess>::iterator itIso;
//   for(itIso=_currentParamJPCIsos12.begin();
//   itIso!=_currentParamJPCIsos12.end(); ++itIso){
//     fitVal.otherParams["Iso12"+(*itIso).first->name()+"Range01"]=(*itIso).second;
//     fitErr.otherParams["Iso12"+(*itIso).first->name()+"Range01"]=0.5;
//    }

//   std::map< std::shared_ptr<const jpcRes>, double,
//   pawian::Collection::SharedPtrLess > currentMagValMap; std::map<
//   std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess >
//   currentPhiValMap; std::map< std::shared_ptr<const jpcRes>, double,
//   pawian::Collection::SharedPtrLess > currentMagErrMap; std::map<
//   std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess >
//   currentPhiErrMap;

//   double magFactor=1./sqrt(_jpcToIGJPCMap.size());
//   for ( std::map< std::shared_ptr<const jpcRes>, std::vector<
//   std::shared_ptr<const IGJPC> >,
//   pawian::Collection::SharedPtrLess>::iterator it = _jpcToIGJPCMap.begin();
//   it!=_jpcToIGJPCMap.end(); ++it){
//     currentMagValMap[it->first] = magFactor;
//     currentMagErrMap[it->first] = magFactor;
//     if(it->first->J > 1./2){
//       currentPhiValMap[it->first] = 0.;
//       currentPhiErrMap[it->first] = 0.3;
//     }
//   }
//   fitVal.MagsJPC["gammap"]=currentMagValMap;
//   fitVal.PhisJPC["gammap"]=currentPhiValMap;
//   fitErr.MagsJPC["gammap"]=currentMagErrMap;
//   fitErr.PhisJPC["gammap"]=currentPhiErrMap;
// }

void gammapBaseLh::fillDefaultParams(
    std::shared_ptr<AbsPawianParameters> fitPar) {
  AbsLh::fillDefaultParams(fitPar);
}

void gammapBaseLh::updateFitParams(
    std::shared_ptr<AbsPawianParameters> fitPar) {
  AbsLh::updateFitParams(fitPar);
}

void gammapBaseLh::fillIsos() {
  // first look for iso0 decay amplitudes and fill it in the map
  std::vector<std::shared_ptr<AbsXdecAmp>>::iterator it;
  for (it = _decAmps.begin(); it != _decAmps.end(); ++it) {
    if (int(2. * (*it)->absDec()->motherIGJPC()->I) == 1) { // I=1/2
      _iso12DecAmps.push_back(*it);
    }
    if (int(2. * (*it)->absDec()->motherIGJPC()->I) == 3) {
      _iso32DecAmps.push_back(*it);
    }
  }

  for (it = _iso12DecAmps.begin(); it != _iso12DecAmps.end(); ++it) {
    _currentParamJPCIsos12[(*it)->jpcPtr()] = 1.;
    _currentParamJPCIsos32[(*it)->jpcPtr()] = 0.;
  }

  for (it = _iso32DecAmps.begin(); it != _iso32DecAmps.end(); ++it) {
    _currentParamJPCIsos32[(*it)->jpcPtr()] = 1.;
    _currentParamJPCIsos12[(*it)->jpcPtr()] = 0.;
  }

  // now look if iso12 and iso32 exisist
  for (it = _iso12DecAmps.begin(); it != _iso12DecAmps.end(); ++it) {
    // find corresponding iso32 partner
    std::vector<std::shared_ptr<AbsXdecAmp>>::iterator itIso1;
    for (itIso1 = _iso32DecAmps.begin(); itIso1 != _iso32DecAmps.end();
         ++itIso1) {
      if ((*(*it)->jpcPtr()) == (*(*itIso1)->jpcPtr())) {
        _currentParamJPCIsos12[(*it)->jpcPtr()] = sqrt(.5);
        _currentParamJPCIsos32[(*it)->jpcPtr()] = sqrt(.5);
      }
    }
  }

  // now look if iso12 and iso32 exisist for the same amplitude
  for (it = _iso12DecAmps.begin(); it != _iso12DecAmps.end(); ++it) {
    std::vector<std::shared_ptr<AbsXdecAmp>>::iterator it2;
    for (it2 = _iso32DecAmps.begin(); it2 != _iso32DecAmps.end(); ++it2) {
      if ((*it)->jpcDecsName() == (*it2)->jpcDecsName()) {
        // found
        _currentParamJPCIsos12[(*it)->jpcPtr()] = sqrt(.7);
        _currentParamJPCIsos32[(*it)->jpcPtr()] = sqrt(.3);
      }
    }
  }
}

// void gammapBaseLh::updateFitParams(fitParCol& theParamVal){
//   AbsLh::updateFitParams(theParamVal);

//   std::map< std::shared_ptr<const jpcRes>, double,
//   pawian::Collection::SharedPtrLess>::iterator itIso;
//   for(itIso=_currentParamJPCIsos12.begin();
//   itIso!=_currentParamJPCIsos12.end(); ++ itIso){
//     double
//     theVal=theParamVal.otherParams.at("Iso12"+(*itIso).first->name()+"Range01");
//     (*itIso).second=theVal;
//     _currentParamJPCIsos32[(*itIso).first]=sqrt(1.-theVal*theVal);
//   }

//   std::map< std::shared_ptr<const jpcRes>, double,
//   pawian::Collection::SharedPtrLess >&
//   magMap=theParamVal.MagsJPC.at("gammap"); std::map< std::shared_ptr<const
//   jpcRes>, double, pawian::Collection::SharedPtrLess >&
//   phiMap=theParamVal.PhisJPC.at("gammap"); std::map< std::shared_ptr<const
//   jpcRes>, std::vector< std::shared_ptr<const IGJPC> >,
//   pawian::Collection::SharedPtrLess>::iterator it; for (it=
//   _jpcToIGJPCMap.begin(); it!=_jpcToIGJPCMap.end(); ++it){
//     double theMag=magMap.at(it->first);
//     _currentParamMags[it->first]=theMag;
//     if(it->first->J > 1./2){
//       double thePhi=phiMap.at(it->first);
//       _currentParamPhis[it->first]=thePhi;
//     }
//   }
// }

void gammapBaseLh::print(std::ostream &os) const {}

void gammapBaseLh::initialize() {
  // check that only one non spin0 final particle exists
  std::vector<Particle *> fsParticles =
      std::static_pointer_cast<GammapChannelEnv>(
          GlobalEnv::instance()->GammapChannel(_channelID))
          ->finalStateParticles();
  std::vector<Particle *>::iterator itParticle;
  bool spinParticleFound = false;

  for (itParticle = fsParticles.begin(); itParticle != fsParticles.end();
       ++itParticle) {
    int current2J = (*itParticle)->twoJ();
    if (current2J > 0.) {
      if (spinParticleFound) {
        Alert << "only one non spin0 final particle allowed!!!" << endmsg;
        exit(0);
      }
      spinParticleFound = true;
    }
  }

  std::shared_ptr<GammapChannelEnv> gammapChannelEnv =
      std::static_pointer_cast<GammapChannelEnv>(
          GlobalEnv::instance()->GammapChannel(_channelID));
  _beamPolarization = GammapBeamPolarization(
      gammapChannelEnv->beamPolFraction(), gammapChannelEnv->beamPolAngle());
  _useReflectivity =
      gammapChannelEnv->parser()->productionFormalism() == "Reflectivity";
  _gammapReactionPtr = gammapChannelEnv->reaction();
  _jpcToIGJPCMap = _gammapReactionPtr->jpcToIGJPCMap();
  _jpcToJPCljMap = _gammapReactionPtr->jpcToJPCLSMap();

  std::vector<std::shared_ptr<AbsDecay>> theDecs =
      gammapChannelEnv->prodDecayList()->getList();
  std::vector<std::shared_ptr<AbsDecay>>::iterator it;
  std::set<std::string> reflectivityProductionChannels;
  for (it = theDecs.begin(); it != theDecs.end(); ++it) {
    std::shared_ptr<AbsXdecAmp> currentAmp;
    if (_useReflectivity) {
      const std::string productionChannel =
          (*it)->daughter1Part()->name() + std::string("_") +
          (*it)->daughter2Part()->name();
      if (!reflectivityProductionChannels.insert(productionChannel).second)
        continue;
      currentAmp = std::shared_ptr<AbsXdecAmp>(
          new ReflectivityDecAmps((*it)->absDecPtr(), _channelID));
    } else {
      currentAmp = XdecAmpRegistry::instance()->getXdecAmp(
          _channelID, (*it)->absDecPtr());
    }
    _decAmps.push_back(currentAmp);
  }
}
