//************************************************************************//
//									  //
//  Copyright 2019 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// TMatrixCompareDynamics class definition file. -*- C++ -*-
// Copyright 2019 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <memory>

#include "PwaUtils/TMatrixDynamics.hh"

class TMatrixCompareDynamics : public TMatrixDynamics{

public:
  TMatrixCompareDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string& pathToConfigParser, std::string& pathToKMatrCompareConfigFile, std::string dataType, std::string projectionParticleNames="");
  virtual ~TMatrixCompareDynamics();

  virtual std::string type() {return "TMatrixDynamics";}
  virtual complex<double> eval(EvtData* theData, AbsXdecAmp* grandmaAmp=0, Spin OrbMom=0);
  virtual void fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillParamNameList();
  virtual bool checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld); 
  virtual void updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar);
  virtual void fillMasses(EvtData* theData);
 
protected:
  int _projectionCompareIndex;
  double _currentOffset;
  std::shared_ptr<TMatrixDynamics> _tMatrDynCompare;
  std::shared_ptr<TMatrixRel> _tMatrCompare;
  std::shared_ptr<KMatrixParser> _kMatrixParserCompare;

  virtual void evalPhaseCompare(EvtData* theData, double currentMass);

private:

};
