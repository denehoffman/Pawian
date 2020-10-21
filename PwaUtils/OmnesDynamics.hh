//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// OmnesDynamics class definition file. -*- C++ -*-
// Copyright 2020 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <memory>

#include "PwaUtils/KMatrixDynamics.hh"

class AbsXdecAmp;
class KMatrixOmnes;
class KPoleOmnes;
class FVectorOmnes;
class KMatrixParser;
class AbsPhaseSpace;
class PVectorBgOmnes;
class AbsPawianParameters;

class OmnesDynamics : public KMatrixDynamics{

public:
  OmnesDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser, ChannelID channelID, std::string projectionParticleNames="");

  virtual ~OmnesDynamics();

  virtual std::string type() {return "OmnesDynamics";}
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom=0);
  
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();
  //  virtual bool checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld); 
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
 
  //  virtual void addGrandMa(std::shared_ptr<AbsDecay> theDec);
  std::string addOneGrandMa(std::string theName);
  //  virtual const std::string& grandMaKey(AbsXdecAmp* grandmaAmp);
  //  virtual const unsigned short grandMaId(AbsXdecAmp* grandmaAmp);
  //virtual std::shared_ptr<FVectorOmnes> fVectorOmnes(std::string);

protected:
  //  std::map<std::string, std::shared_ptr<PVectorBgOmnes> > _pVecMap;

  std::shared_ptr<PVectorBgOmnes> makeNewPVecOmnes();
  virtual void init();

  std::shared_ptr<KMatrixOmnes> _kMatrOmnes;
  std::vector< std::shared_ptr<KPoleOmnes> > _kPolesOmnes;
private:

};
