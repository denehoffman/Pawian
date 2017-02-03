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

#include "PwaUtils/TMatrixDynamics.hh"

class AbsXdecAmp;
class KMatrixRel;
class KPole;
class FVector;
class PVectorKPiSFocus;
class KMatrixParser;
class AbsPhaseSpace;
class PVectorRel;
class AbsPawianParameters;

class KMatrixDynamics : public TMatrixDynamics{

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
  int _orderPVecBg;
  double _currentMass;
  std::map<std::string, std::shared_ptr<AbsXdecAmp> >  _grandMaAmpMap;
  
  std::map<std::string, std::shared_ptr<FVector> > _fVecMap;
  std::map<std::string, std::shared_ptr<PVectorRel> > _pVecMap;

  std::map<std::string, std::map<std::string, double> > _currentbFactorMap;
  std::vector< std::vector<double> >  _currentPVecBgTerms;
  std::vector< std::vector<std::string> >  _bgPVecTermNames;
  
  CacheVector<std::map<std::string, complex<float> > > _cachedStringMap;
  CacheVector<std::map<std::string, std::map<Spin, complex<float> > > > _cachedStringOrbMap;
  std::map<int, std::map<std::string, bool > > _alreadyCached;

  virtual void init();
  std::shared_ptr<PVectorRel> makeNewPVec();
private:

};
