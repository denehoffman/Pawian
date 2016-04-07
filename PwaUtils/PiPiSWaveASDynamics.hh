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

// PiPiSWaveASDynamics class definition file. -*- C++ -*-
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
class FVector;
class ParticleTable;
class KMatrixPiPiS;
class PVectorSlowCorRel;
class AbsPhaseSpace;
class KPole;
class AbsPawianParameters;

class PiPiSWaveASDynamics : public AbsDynamics{

public:
  PiPiSWaveASDynamics(std::string& massKey, std::vector<Particle*>& fsParticles, Particle* mother, ParticleTable* thePdtTable);
  virtual ~PiPiSWaveASDynamics();

  virtual std::string type() {return "PiPiSWaveASDynamics";}
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom=0);
  
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList(); 
  virtual bool checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld);
  
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar); 
  virtual void addGrandMa(std::shared_ptr<AbsDecay> theDec);
  virtual const std::string& grandMaKey(AbsXdecAmp* grandmaAmp);
  virtual const unsigned short grandMaId(AbsXdecAmp* grandmaAmp);

protected:
  virtual int projectionIndex(std::vector<Particle*>& fsParticles);

  double _currentMass;
  std::map<std::string, std::shared_ptr<AbsXdecAmp> >  _grandMaAmpMap;

  std::map<std::string, std::shared_ptr<FVector> > _fVecMap;

  std::map<std::string, std::map<std::string, double> > _currentbFactorMagMap;
  std::map<std::string, std::map<std::string, double> > _currentbFactorPhiMap;
  std::map<std::string, std::map<std::string, double> > _currentfProdMagMap;
  std::map<std::string, std::map<std::string, double> > _currentfProdPhiMap;
  std::map<std::string, double > _currentS0Map; 

  CacheVector<std::map<std::string, complex<float> > > _cachedStringMap;
  std::map<int, std::map<std::string, bool > > _alreadyCached;
  std::map<std::string, bool > _recalcMap;

  ParticleTable* _pdtTable;
  int _projectionIndex;
  std::shared_ptr<KMatrixPiPiS> _KmatrixPiPiS;
  std::map<std::string, std::shared_ptr<PVectorSlowCorRel> > _pVecMap;
  vector<std::shared_ptr<AbsPhaseSpace> > _phpVec;
  vector<std::shared_ptr<KPole> > _kPoles;
  std::map<std::string, std::vector<std::string> > _paramNameListMap;
  
private:
  std::shared_ptr<PVectorSlowCorRel> makeNewPVec();
};
