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

// ProdChannelInfo class definition file. -*- C++ -*-
// Copyright 2015 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>

class Particle;

class ProdChannelInfo {

public:
  ProdChannelInfo();
  ProdChannelInfo(std::string& stringFromParser); //constructor by making use of the string provided by the parser (very specific!!!)
  virtual ~ProdChannelInfo();
  const bool isProductionChannel() const {return _isProdChannel;}
  const bool isFormation() const {return _isFormation;}
  std::pair<Particle*, Particle*> productionPair() const {return _prodPair;}
  Particle*  formationParticle() const {return _formationParticle;}
  const bool withProdBarrier() const {return _withProBarrier;}
  const std::string prodBarrierType() const {return _prodBarrierType;}
  const double qRPod() const {return _qRProd;}
  const double m0decRadM1Prod() const {return _m0decRadM1Production;}
  const std::string formationDynType() const {return _formationDynType;}
  const std::string additionalFormationDynInfos() const {return _additionalFormationDynInfos;}
  virtual void print(std::ostream& os) const;

 protected:
  bool _isProdChannel;
  bool _isFormation;
  std::pair<Particle*, Particle*> _prodPair;
  Particle* _formationParticle;
  bool _withProBarrier;
  std::string _prodBarrierType;
  std::string _formationDynType;
  std::string _additionalFormationDynInfos;
  double _qRProd;
  double _m0decRadM1Production;
};
