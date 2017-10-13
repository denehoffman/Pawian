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

// AbsDecay class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <memory>

#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "Utils/PawianCollectionUtils.hh"

class Particle;
class EvtData;
class AbsDynamics;
class ProdChannelInfo;

class AbsDecay : public std::enable_shared_from_this<AbsDecay>{

public:
  typedef enum decayLevel {noLevel=-1, isProdAmp, firstLevel, secondLevel, thirdLevel, fourthLevel, fifthLevel, severalLevels} decLevel;
  
  AbsDecay(Particle* mother, Particle* daughter1, Particle* daughter2, ChannelID channelId, std::string typeName="AbsDecay");
  AbsDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, Particle* daughter2, std::string motherName, ChannelID channelId, std::string typeName="AbsDecay");
  AbsDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, std::string motherName, ChannelID channelId, std::string typeName="AbsDecay");

  virtual ~AbsDecay();

  //  virtual AbsDecay* clone_() const = 0;
  virtual const std::string name() const {return _name;}
  const unsigned short nameId() const {return _nameId;}
  const std::string wignerDKey() {return _wignerDKey;}
  const unsigned short wignerDId() {return _wignerDId;}
  const std::string wignerDRefKey() {return _wignerDRefKey;}
  const unsigned short wigDWigDRefId() {return _wigDWigDRefId;}
  const std::string refKey() {return _refKey;}
  virtual std::string fitParSuffix() const {return _fitParamSuffix;}
  void setFitParSuffix(std::string& suffix) {_fitParamSuffix = suffix;}
  virtual std::string& massParKey() {return _massParamKey;}
  virtual unsigned short massParId() {return _massParamId;}
  void setMassParKey(const std::string& newKey);
  virtual std::string& prodParKey() {return _prodParamKey;}
  void setProdParKey(const std::string& newKey) {_prodParamKey = newKey;}
  void setPreFactor(double thePreFactor) {_preFactor=thePreFactor;}
  double preFactor() {return _preFactor;}
  std::shared_ptr<const jpcRes> motherJPC(){ return _motherJPCPtr;}
  std::shared_ptr<const IGJPC> motherIGJPC(){ return _motherIGJPCPtr;}
  //  std::vector< std::shared_ptr<const JPCLS> > JPCLSAmps(){ return _JPCLSDecAmps;}
  std::shared_ptr<AbsDecay> decDaughter1() {return _absDecDaughter1;}
  std::shared_ptr<AbsDecay> decDaughter2() {return _absDecDaughter2;}
  bool hasMother() {return _hasMotherPart;}
  bool isDaughter1Stable() {return _daughter1IsStable;}
  bool isDaughter2Stable() {return _daughter2IsStable;}
  virtual bool isTensorAmp() {return false;}
  virtual void setDecayLevel(decLevel theLevel);
  virtual void setDecayLevelTree(decLevel theLevel, std::shared_ptr<AbsDecay> motherDecPtr, std::shared_ptr<AbsDecay> prodDecPtr); 
  virtual decLevel whichDecayLevel() {return _decLevel;}

  std::vector<Particle*> finalStateParticles() {return _finalStateParticles;}
  std::vector<Particle*> finalStateParticlesDaughter1() {return _finalStateParticlesDaughter1;}
  std::vector<Particle*> finalStateParticlesDaughter2() {return _finalStateParticlesDaughter2;}
  std::vector< std::shared_ptr<AbsDecay> > prodAmpRefList() {return _prodAmpRefList;}
  std::vector< std::shared_ptr<AbsDecay> > motherAmpRefList() {return _motherAmpRefList;}
  std::vector<std::string> fsParticleNameList() {return _fsParticleNameList;}
  virtual void addToFsParticleNameList(const std::string& name);

  //  virtual void fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, Vector4<double>& prodParticle4Vec, EvtData* evtData, std::string& refKey);
  virtual void fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, Vector4<double>& prodParticle4Vec, EvtData* evtData);
  void enableDynamics(std::string& dynString, std::vector<std::string>& additionalStringVec);
  std::shared_ptr<AbsDynamics> getDynamics(){return _absDynPtr;}
  double barrierqR() {return _qR;}
  virtual void print(std::ostream& os) const;

  Particle* motherPart() {return _mother;}
  Particle* daughter1Part() {return _daughter1;}
  Particle* daughter2Part() {return _daughter2;}

  virtual std::string dynType() {return _dynType;}
  //  virtual std::string& dynKey() {return _dynKey;}
  virtual std::pair<Particle*, Particle*>& firstDecayChannel() { return _decPair1stChannel;}
  virtual std::pair<Particle*, Particle*>& secondDecayChannel() { return _decPair2ndChannel;}
  virtual std::shared_ptr<AbsDecay> absDecPtr() {return shared_from_this();}
  virtual std::string type() {return _typeName;}
  virtual void extractStates() =0;

  double isospinCG() {return _isospinClebschG;}

  Spin iDaughter1() {return _idaughter1;}
  Spin i3Daughter1() {return _i3daughter1;}
  Spin iDaughter2() {return _idaughter2;}
  Spin i3Daughter2() {return _i3daughter2;}
  virtual void disableIsospin(){_useIsospin=false;}
  bool useIsospin(){ return _useIsospin;}
  bool isWeakDeacy() {return _isWeakDecay;}
  virtual std::string pathToConfigParser() {return _pathParserFile;}
  virtual std::string projectionParticleNames() {return _projectionParticleNames;}
  std::shared_ptr<ProdChannelInfo> prodChannelInfo() {return _prodChannelInfo;}
  void setProdChannelInfo(std::shared_ptr<ProdChannelInfo> prodChannelInfo) { _prodChannelInfo=prodChannelInfo;}
  void setProductionAmp() {_isProdAmp=true;}
  bool isProductionAmp() {return _isProdAmp;}
  virtual void setWigDRefKey(std::string& ref);

  void resetFilledMap();
  virtual void enableProdBarrier();
  bool useProdBarrier() {return _useProdBarrier;}
  double massSumFsParticles() {return _massSumFsParticles;}
  double massSumFsParticlesDec1();
  double massSumFsParticlesDec2();
  int orbMomMin() {return _Lmin;}
  ChannelID channelId() {return _channelId;}
protected:
  virtual Vector4<double> beamVecCollider(Vector4<double>& sqrts, double massBeam); 
  virtual void extractLmin();
  std::string _typeName;
  ChannelID _channelId;

  Particle* _mother;
  Particle* _daughter1;
  Particle* _daughter2;

  bool _daughter1IsStable;
  bool _daughter2IsStable;
  bool _hasMotherPart;

  std::shared_ptr<const jpcRes> _motherJPCPtr;
  std::shared_ptr<const IGJPC> _motherIGJPCPtr;
  std::shared_ptr<const IGJPC> _daughter1IGJPCPtr;
  std::shared_ptr<const IGJPC> _daughter2IGJPCPtr;

  double _isospinClebschG;
  double _qR;

  std::string _name;
  unsigned short _nameId;
  std::string _fitParamSuffix;
  std::string _massParamKey;
  unsigned short _massParamId;
  std::string _prodParamKey;
  std::string _wignerDKey;
  unsigned short _wignerDId;
  unsigned short _wignerDqId;
  unsigned short _wignerDqNormId;
  std::string _refKey;
  std::string _wignerDRefKey;
  unsigned short _wigDWigDRefId;
  std::string _dynType;
  bool _dynEnabled;

  double _preFactor;
  std::string _pathParserFile;
  std::string _projectionParticleNames;

  std::vector<std::string> _fsParticleNameList;
  std::vector< std::shared_ptr<AbsDecay> > _prodAmpRefList;
  std::vector< std::shared_ptr<AbsDecay> > _motherAmpRefList;
  std::shared_ptr<AbsDecay> _absDecDaughter1;
  std::shared_ptr<AbsDecay> _absDecDaughter2;

  std::vector<Particle*> _finalStateParticles;
  std::vector<Particle*> _finalStateParticlesDaughter1;
  std::vector<Particle*> _finalStateParticlesDaughter2;
  std::pair<Particle*, Particle*> _decPair1stChannel;
  std::pair<Particle*, Particle*> _decPair2ndChannel;

  std::shared_ptr<AbsDynamics> _absDynPtr;
  std::map<int, bool> _alreadyFilledMap;

  Spin _idaughter1;
  Spin _i3daughter1;
  Spin _idaughter2;
  Spin _i3daughter2;

  bool _isDaughter1Photon;
  bool _isDaughter2Photon;

  int _gParity;

  bool _useIsospin;
  bool _isWeakDecay;
  bool _isProdAmp;
  bool _useProdBarrier;
  double _massSumFsParticles;

  int _Lmin;
  decLevel _decLevel;
  std::string _motherName;
  std::shared_ptr<ProdChannelInfo> _prodChannelInfo;
};
