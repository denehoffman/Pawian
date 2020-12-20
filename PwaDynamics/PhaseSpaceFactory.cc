//************************************************************************//
//									  //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// PhaseSpaceFactory class definition file. -*- C++ -*-
// Copyright 2016 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaDynamics/PhaseSpaceFactory.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceIsobar.hh"
#include "PwaDynamics/PhaseSpaceIsobarReid.hh"
#include "PwaDynamics/PhaseSpaceIsobarReidAngularMomentum.hh"
#include "PwaDynamics/PhaseSpaceIsobarDudek.hh"
#include "PwaDynamics/PhaseSpaceIsobarDudekUnstableRhoPi.hh"
//#include "PwaDynamics/PhaseSpaceIsobarBBUnstable.hh"
#include "PwaDynamics/PhaseSpaceIsobarLUT.hh"
#include "PwaDynamics/PhaseSpaceIsobarAS.hh"
#include "PwaDynamics/PhaseSpace4Pi.hh"
#include "ErrLogger/ErrLogger.hh"


PhaseSpaceFactory* PhaseSpaceFactory::_instance=0;

PhaseSpaceFactory* PhaseSpaceFactory::instance()
{
  if (0==_instance) _instance = new PhaseSpaceFactory();
  return _instance;
}

PhaseSpaceFactory::PhaseSpaceFactory()
{
}

PhaseSpaceFactory::~PhaseSpaceFactory()
{
}

std::shared_ptr<AbsPhaseSpace> PhaseSpaceFactory::getPhpPointer(std::string type, std::vector<double> masses){
  std::shared_ptr<AbsPhaseSpace> result;
  //masses should constists of at least two entries
  if(masses.size() <2){
   Alert << "number of final state particles must be at least 2!!!" << endmsg;
   exit(1);
  }
  
  double mass1 = masses.at(0);
  double mass2 = masses.at(1);

  if(type=="Default"){
	 result= std::shared_ptr<AbsPhaseSpace>(new PhaseSpaceIsobar(mass1, mass2));
 }
 else if(type=="Reid"){
	 result= std::shared_ptr<AbsPhaseSpace>(new PhaseSpaceIsobarReid(mass1, mass2));
 }
 else if(type=="ReidAngularMomentum"){
   result= std::shared_ptr<AbsPhaseSpace>(new PhaseSpaceIsobarReidAngularMomentum(mass1, mass2));
 }
 else if(type=="Dudek"){
	 result= std::shared_ptr<AbsPhaseSpace>(new PhaseSpaceIsobarDudek(mass1, mass2));
 }
 else if(type=="DudekUnstableRhoPi"){
   result= std::shared_ptr<AbsPhaseSpace>(new PhaseSpaceIsobarDudekUnstableRhoPi(mass1, mass2));
 }
 else if(type=="AS"){
	 result= std::shared_ptr<AbsPhaseSpace>(new PhaseSpaceIsobarAS(mass1, mass2));
 }
 else if(type=="4pi"){
	 result= std::shared_ptr<AbsPhaseSpace>(new PhaseSpace4Pi());
 }
// else if(type.substr (0,10)=="BBUnstable"){
//	 result= std::shared_ptr<AbsPhaseSpace>(new PhaseSpaceIsobarBBUnstable(mass1, mass2, type));
// }
 else if(type.substr (0,3)=="LUT"){
	 result= std::shared_ptr<AbsPhaseSpace>(new PhaseSpaceIsobarLUT(mass1, mass2, type));
 }
 else{
    Alert << "phase space description of the type:\t" << type << "\tdoes not exist" << endmsg;
    exit(1);
  }

  return result;
}

