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

// GammapChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf, Julian Pychy

#include "gammapUtils/GammapChannelEnv.hh"
#include "ConfigParser/gammapParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "gammapUtils/gammapHist.hh"
#include "gammapUtils/gammapReaction.hh"

GammapChannelEnv::GammapChannelEnv(gammapParser *theParser)
    : AbsChannelEnv(theParser, AbsChannelEnv::CHANNEL_GAMMAP), _lmax(0),
      _beamPolFraction(0.), _beamPolAngle(0.), _theGamPParser(theParser) {}

void GammapChannelEnv::setupChannel(ChannelID id) {

  AbsChannelEnv::setupGlobal(id);

  double pMass =
      GlobalEnv::instance()->particleTable()->particle("proton")->mass();

  // Use a nominal GlueX coherent-peak photon energy if generation has no
  // configured total initial-state four-vector.
  if (_initial4Vec == Vector4<double>(0., 0., 0., 0.)) {
    double nominalPhotonEnergy = 7.;
    _initial4Vec = Vector4<double>(pMass + nominalPhotonEnergy, 0., 0.,
                                   nominalPhotonEnergy);
    WarningMsg << "NO initial 4-vector set in config file! Using px="
               << _initial4Vec.Px() << ", py=" << _initial4Vec.Py()
               << ", pz=" << _initial4Vec.Pz() << ", E=" << _initial4Vec.E()
               << " instead!" << endmsg;
  }
  _cmEnergy = _initial4Vec.M();

  // Lmax
  _lmax = _theGamPParser->getLMax();
  _beamPolFraction = _theGamPParser->beamPolFraction();
  _beamPolAngle = _theGamPParser->beamPolAngle();
  if (_beamPolFraction < 0. || _beamPolFraction > 1.) {
    Alert << "beamPolFraction must be in [0, 1], got " << _beamPolFraction
          << endmsg;
    exit(1);
  }

  // individual Lmax settings
  std::vector<std::string> theDropGammapLForParticles =
      _theGamPParser->dropGammapLForParticle();
  for (auto ldropIt = theDropGammapLForParticles.begin();
       ldropIt != theDropGammapLForParticles.end(); ++ldropIt) {
    std::string particle;
    short l;
    std::stringstream stringStr;
    stringStr << (*ldropIt);
    stringStr >> particle >> l;
    _dropGammapLForParticleData[particle].push_back(l);
  }

  // gammap reaction
  _gammapReaction = std::shared_ptr<gammapReaction>(
      new gammapReaction(_prodChannelInfoList, id, _lmax));

  // fill prodDecayList
  std::vector<std::string> additionalStringVecDummy;
  std::string dynTypeDefault = "WoDynamics";

  if (_theGamPParser->productionFormalism() == "Cano" ||
      _theGamPParser->productionFormalism() == "Reflectivity") {
    std::vector<std::shared_ptr<IsobarLSDecay>> prodDecs =
        _gammapReaction->productionDecays();
    std::vector<std::shared_ptr<IsobarLSDecay>>::iterator itDec;
    for (itDec = prodDecs.begin(); itDec != prodDecs.end(); ++itDec) {
      if (_theGamPParser->productionFormalism() == "Cano" &&
          (*itDec)->prodChannelInfo()->withProdBarrier())
        (*itDec)->enableProdBarrier();
      else
        (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      _prodDecList->addDecay(*itDec);
    }
  } else if (_theGamPParser->productionFormalism() == "Tensor") {
    std::vector<std::shared_ptr<IsobarTensorDecay>> prodDecs =
        _gammapReaction->productionTensorDecays();
    std::vector<std::shared_ptr<IsobarTensorDecay>>::iterator itDec;
    for (itDec = prodDecs.begin(); itDec != prodDecs.end(); ++itDec) {
      (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      _prodDecList->addDecay(*itDec);
    }
  } else if (_theGamPParser->productionFormalism() == "Heli") {
    std::vector<std::shared_ptr<IsobarHeliDecay>> prodDecs =
        _gammapReaction->productionHeliDecays();
    std::vector<std::shared_ptr<IsobarHeliDecay>>::iterator itDec;
    for (itDec = prodDecs.begin(); itDec != prodDecs.end(); ++itDec) {
      if ((*itDec)->prodChannelInfo()->withProdBarrier())
        (*itDec)->enableProdBarrier();
      else
        (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      _prodDecList->addDecay(*itDec);
    }
  } else {
    Alert << "production formalism\t" << _theGamPParser->productionFormalism()
          << "\t is not supported!!!" << endmsg;
    exit(0);
  }

  AbsChannelEnv::setPrefactors();
  AbsChannelEnv::replaceParameterSuffixes();
  AbsChannelEnv::replaceMassKeys();
  AbsChannelEnv::addDynamics();
  AbsChannelEnv::setDecayLevels();

  // spin density particles
  _spinDensity = _theGamPParser->spinDensityNames();
}

std::shared_ptr<AbsHist>
GammapChannelEnv::CreateHistInstance(std::string additionalSuffix) {

  return std::shared_ptr<AbsHist>(new gammapHist(additionalSuffix));
}
