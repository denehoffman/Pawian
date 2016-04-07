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

// KMatrixDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

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
class KPole;
class FVector;
class PVectorKPiSFocus;
class KMatrixParser;
class AbsPhaseSpace;
class PVectorRel;
class AbsPawianParameters;

class KMatrixDynamics : public AbsDynamics{

public:
  KMatrixDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser);
  virtual ~KMatrixDynamics();

  virtual std::string type() {return "KMatrixDynamics";}
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom=0);
  
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();
  virtual bool checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld); 
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
 
  virtual void addGrandMa(std::shared_ptr<AbsDecay> theDec);
  virtual const std::string& grandMaKey(AbsXdecAmp* grandmaAmp);
  virtual const unsigned short grandMaId(AbsXdecAmp* grandmaAmp);

protected:
  std::string _kMatName;
  int _projectionIndex;
  int _orderKMatBg;
  bool _withKMatAdler;
  double _currentMass;
  std::map<std::string, std::shared_ptr<AbsXdecAmp> >  _grandMaAmpMap;
  
  std::shared_ptr<KMatrixRel> _kMatr;
  std::vector< std::shared_ptr<KPole> > _kPoles;
  std::map<std::string, std::shared_ptr<FVector> > _fVecMap;
  std::map<std::string, std::shared_ptr<PVectorRel> > _pVecMap;
  std::vector<std::shared_ptr<AbsPhaseSpace> > _phpVecs;
  std::vector< std::string> _poleNames;
  std::vector< std::string> _gFactorNames;

  std::map<std::string, std::map<std::string, double> > _currentbFactorMap;
  std::vector<double> _currentPoleMasses;
  std::map<int, std::vector<double> > _currentgFactorMap;
  std::vector< std::vector< std::vector<double> > > _currentBgTerms;
  std::vector< std::vector< std::vector<std::string> > > _bgTermNames;
  double _currentAdler0;
  
  CacheVector<std::map<std::string, complex<float> > > _cachedStringMap;
  CacheVector<std::map<std::string, std::map<Spin, complex<float> > > > _cachedStringOrbMap;
  std::map<int, std::map<std::string, bool > > _alreadyCached;
  std::map<std::string, bool > _recalcMap;
  std::shared_ptr<KMatrixParser> _kMatrixParser;
  std::map<std::string, std::vector<std::string> > _paramNameListMap;

  virtual void init();
  std::shared_ptr<PVectorRel> makeNewPVec();
private:

};
