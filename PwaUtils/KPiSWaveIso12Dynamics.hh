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

// KPiSWaveIso12Dynamics class definition file. -*- C++ -*-
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
class KMatrixKPiSFocus;
class FVector;
class PVectorKPiSFocus;

class KPiSWaveIso12Dynamics : public AbsDynamics{

public:
  KPiSWaveIso12Dynamics(std::string& massKey, std::vector<Particle*>& fsParticles, Particle* mother);
  virtual ~KPiSWaveIso12Dynamics();

  virtual std::string type() {return "KPiSWaveIso12Dynamics";}
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom=0);
  
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual bool checkRecalculation(fitParams& theParamVal);
  virtual void updateFitParams(fitParams& theParamVal);
  virtual void addGrandMa(std::shared_ptr<AbsDecay> theDec);
  virtual const std::string& grandMaKey(AbsXdecAmp* grandmaAmp);

protected:
  //  std::string _massKey;
  double _currentMass;
  std::map<std::string, std::shared_ptr<AbsXdecAmp> >  _grandMaAmpMap;

  std::shared_ptr<KMatrixKPiSFocus> _kMatr;
  std::map<std::string, std::shared_ptr<FVector> > _fVecMap;
  std::map<std::string, std::shared_ptr<PVectorKPiSFocus> > _pVecMap;

  std::map<std::string, std::map<std::string, double> > _currentbFactorMap;
  std::map<std::string, std::map<std::string, double> > _currentaProdMap; 
  std::map<std::string, std::map<std::string, double> > _currentbProdMap; 
  std::map<std::string, std::map<std::string, double> > _currentcProdMap;
  std::map<std::string, std::map<std::string, double> > _currentphaseProdMap;

  std::map<int, std::map<std::string, complex<float> > > _cachedStringMap;
  std::map<std::string, bool > _recalcMap;
private:

};
