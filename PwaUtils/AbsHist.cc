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

// AbsHist class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "PwaUtils/AbsHist.hh"

#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/IdStringMapRegistry.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "FitParams/AbsPawianParameters.hh"

#include "ErrLogger/ErrLogger.hh"

AbsHist::AbsHist(std::string additionalSuffix, bool withTruth) :
  _fillTruths(withTruth)
  ,_additionalSuffix(additionalSuffix)
  ,_weightToWrite(1.)
  ,_fsParticles(GlobalEnv::instance()->Channel()->finalStateParticles())
  ,_cmMass(GlobalEnv::instance()->Channel()->cmEnergy())
{
  _angleHistDataVec=GlobalEnv::instance()->Channel()->angleHistDataVec();
  _angleHistDataVec2D=GlobalEnv::instance()->Channel()->angleHistDataVec2D();
}

AbsHist::~AbsHist(){
}

