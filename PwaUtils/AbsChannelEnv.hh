//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

// AbsChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#pragma once

#include <vector>
#include <memory>
#include "PwaUtils/AbsHist.hh"
#include "PwaUtils/WaveContribution.hh"
#include "PwaUtils/DataUtils.hh"

class Particle;
class AbsDecayList;
class IsobarDecayList;
class ParserBase;
class AbsLh;
class MassRangeCut;
class PhpGenDynamics;
class ProdChannelInfo;
class FsParticleProjections;

typedef short ChannelID;


class AbsChannelEnv
{
public:

  //  virtual void setup(ChannelID id);
  virtual void setupChannel(ChannelID id)=0;
  const int noFinalStateParticles() {return _noFinalStateParticles;}
  std::vector<Particle*> finalStateParticles() {return _finalStateParticles;}
  std::shared_ptr<AbsDecayList> absDecayList() {return _absDecList;}
  std::shared_ptr<AbsDecayList> prodDecayList() {return _prodDecList;}
  bool checkReactionChain();
  void setWignerDRefs();
  std::vector< std::shared_ptr<ProdChannelInfo> > prodChannelInfoList() {return _prodChannelInfoList;}
  std::vector<std::vector<std::string> >& histMassSystems() {return _histMassSystems;}
  std::vector<std::shared_ptr<angleHistData> >& angleHistDataVec() {return _angleHistDataVec;}
  std::vector<std::shared_ptr<angleHistData2D> >& angleHistDataVec2D() {return _angleHistDataVec2D;}
  std::vector<std::shared_ptr<calcContributionData> >& calcContributionDataVec() {return _calcContributionDataVec;}
  const bool useMassRange() const {return _useMassRange;}
  const bool usePhpGenDynamics() const {return _usePhpGenDynamics;}
  bool useDataEvtWeight(){return _useDataEvtWeight;}
  bool useMCEvtWeight(){return _useMCEvtWeight;}
  std::vector< std::shared_ptr<MassRangeCut> > massRangeCuts() {return _massRangeCuts;}
  std::vector< std::shared_ptr<PhpGenDynamics> > phpGenDynamics() {return _phpGenDynamics;}
  Vector4<double>&  initial4Vec() {return _initial4Vec;}
  ParserBase* parser(){ return _theParser;}
  void SetLh(std::shared_ptr<AbsLh> theLh){ _theLh = theLh; }
  std::shared_ptr<AbsLh> Lh();
  ChannelID channelID(){return _channelID;}
  short channelType(){return _channelType;}
  virtual const std::string  channelTypeName()=0;
  virtual std::shared_ptr<AbsHist> CreateHistInstance(std::string additionalSuffix="", bool withTruth=false);
  std::shared_ptr<AbsPawianParameters> defaultPawianParams();
  std::shared_ptr<FsParticleProjections> getFsParticleProjectionsPtr() {return _fsParticleProjections;}
  double cmEnergy() {return _cmEnergy;}

  static short CHANNEL_PBARP;
  static short CHANNEL_EPEM;
  static short CHANNEL_RES;
  static short CHANNEL_GAMMAP;
  static short CHANNEL_PIPISCATTERING;
  
protected:
   int _channelID;
   short _channelType;
   bool _alreadySetUp;
   bool _useDataEvtWeight;
  bool _useMCEvtWeight;
  AbsChannelEnv(ParserBase* theParser, short channelType);
  virtual void setupGlobal(ChannelID id);
  virtual void replaceParameterSuffixes();
  virtual void replaceMassKeys();
  virtual void addDynamics();
  virtual void setDecayLevels();
  virtual void setPrefactors();
   int _noFinalStateParticles;
   std::vector<Particle*> _finalStateParticles;
   std::shared_ptr<AbsDecayList> _absDecList;
   std::shared_ptr<AbsDecayList> _prodDecList;
  std::vector< std::shared_ptr<ProdChannelInfo> > _prodChannelInfoList;
  std::shared_ptr<FsParticleProjections> _fsParticleProjections;
   std::vector<std::vector<std::string> > _histMassSystems;
   std::vector<std::shared_ptr<angleHistData> > _angleHistDataVec;
   std::vector<std::shared_ptr<angleHistData2D> > _angleHistDataVec2D;
   std::vector<std::shared_ptr<calcContributionData> > _calcContributionDataVec;
   bool _useMassRange;
  bool _usePhpGenDynamics;
   double _massMin;
   double _massMax;
   Vector4<double> _initial4Vec;
   std::map<std::string, double> _preFactorMap;
   ParserBase* _theParser;
   std::shared_ptr<AbsLh> _theLh;
  std::vector< std::shared_ptr<MassRangeCut> > _massRangeCuts;
  std::vector< std::shared_ptr<PhpGenDynamics> > _phpGenDynamics;
  double _cmEnergy;
};
