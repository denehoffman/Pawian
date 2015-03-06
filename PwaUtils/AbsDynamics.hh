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

// AbsDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <mutex>
#include <memory>

#include "PwaUtils/EvtDataBaseList.hh"
#include "FitParams/FitParColBase.hh"
#include "PwaUtils/AbsParamHandler.hh"

class Particle;
class AbsDecay;
class AbsXdecAmp;

class AbsDynamics : public AbsParamHandler{

public:
  AbsDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother );
  virtual ~AbsDynamics();

  virtual std::string name() {return _name;};
  virtual std::string type()=0;
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom=0)=0;
  virtual void cacheAmplitudes();
  virtual void fillMasses(EvtData* theData);
  virtual void setMassKey(std::string& theMassKey){_massKey=theMassKey;}
  virtual void addGrandMa(std::shared_ptr<AbsDecay> theDec) {;}
  virtual const std::string& grandMaKey(AbsXdecAmp* grandmaAmp){return _grandmaKey;}
  virtual bool isLdependent() {return _isLdependent;}

protected:
  std::string _name;
  std::string _massKey;
  std::vector<Particle*> _fsParticles;
  Particle* _mother;
  CacheVector<complex<float> >  _cachedMap;
  std::mutex theMutex;
  std::string _dynKey;
  std::string _grandmaKey;
  bool _isLdependent;
private:

};
