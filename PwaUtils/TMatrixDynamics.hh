//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// TMatrixDynamics class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <memory>

#include "PwaUtils/AbsDynamics.hh"

class AbsXdecAmp;
class KMatrixRel;
class TMatrixRel;
class KPole;
class KMatrixParser;
class AbsPhaseSpace;
class AbsPawianParameters;

class TMatrixDynamics : public AbsDynamics{

public:
  TMatrixDynamics(std::shared_ptr<KMatrixParser> kMatrixParser);
  //  TMatrixDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser, std::string dataType="Elasticity");
  TMatrixDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser, std::string dataType, std::string projectionParticleNames="");
  TMatrixDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother);
  virtual ~TMatrixDynamics();

  virtual std::string type() {return "TMatrixDynamics";}
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp=0, Spin OrbMom=0);
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();
  virtual bool checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld); 
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);

  virtual std::shared_ptr<KMatrixRel> getKMatix() {return _kMatr;}
  virtual std::shared_ptr<TMatrixRel> getTMatix() {return _tMatr;}
  virtual std::vector< std::string>   gFactorNames() {return _gFactorNames;}
  virtual std::string KMatixName() {return _kMatName;}
  virtual int orbitalL() {return _orbitalL;} 
  virtual void fillMasses(EvtData* theData);
  virtual int decProjectionIndex() {return _decProjectionIndex;}
  int prodProjectionIndex() {return _prodProjectionIndex;}
  void setProdProjectionIndex(int idx);
  void doTcheck(Spin OrbMom);
  
protected:
  std::string _kMatName;
  int _prodProjectionIndex;
  int _decProjectionIndex;
  int _orderKMatBg;
  bool _withKMatAdler;
  short _dataTypeID;
  bool _prodIsNotDecChannel;
  std::string _paramNameRelPhase;
  double _currentRelPhase;
  int _orbitalL; 
  
  std::shared_ptr<KMatrixRel> _kMatr;
  std::shared_ptr<TMatrixRel> _tMatr;
  std::vector< std::shared_ptr<KPole> > _kPoles;
  std::vector<std::shared_ptr<AbsPhaseSpace> > _phpVecs;
  std::vector< std::string> _poleNames;
  std::vector< std::string> _gFactorNames;

  std::vector<double> _currentPoleMasses;
  std::map<int, std::vector<double> > _currentgFactorMap;
  std::vector< std::vector< std::vector<double> > > _currentBgTerms;
  std::vector< std::vector< std::vector<std::string> > > _bgTermNames;
  double _currentAdler0;

  std::map<std::string, bool > _recalcMap;  
  std::shared_ptr<KMatrixParser> _kMatrixParser;
  std::map<std::string, std::vector<std::string> > _paramNameListMap;
  std::map<unsigned int, double > _noRotationMap;
  std::string _projectionParticleNames;
  std::string _pathToKMatCfgParser;
  
  virtual void init();
  virtual void evalElasticity(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalPhase(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalRelativePhase(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalArgandUnits(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalTreal(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalTimag(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalTabs(EvtData* theData, double currentMass, Spin OrbMom);
  virtual void evalPhasePhi(EvtData* theData, double currentMass, Spin OrbMom);
  unsigned int noOfRotations(double currentMass);
  
private:

};
