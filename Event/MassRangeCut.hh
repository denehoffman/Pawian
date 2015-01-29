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

#pragma once

#include "ErrLogger/ErrLogger.hh"
#include <string>
#include <vector>

class Particle;

class MassRangeCut
{
public:
  MassRangeCut(std::string& rangeAndParticleNames, std::vector<Particle*>& finalStateParticles);
  ~MassRangeCut();

  const double massMin() const {return _massMin;}
  const double massMax() const {return _massMax;}
  bool isMassRangeParticle(Particle* compParticle);
  std::vector<unsigned int> particleIds() {return _iDparticlesMassRange;}

private:
  double _massMin;
  double _massMax;
  std::vector<unsigned int> _iDparticlesMassRange;
  std::vector<Particle*> _particlesMassRange;
};
