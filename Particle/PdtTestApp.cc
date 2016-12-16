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

#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Particle/PdtParser.hh"

#include "ErrLogger/ErrLogger.hh"

#include <iostream>
#include <string>
#include <cstdlib>

int main()
{

  ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);

  ParticleTable pTable;
  PdtParser parser;
  std::string theSourcePath=getenv("TOP_DIR"); 
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");

  if (parser.parse(pdtFile, pTable)) {
    pTable.print(std::cout);

    std::string name("omega");
    Particle* omega = pTable.particle(name);
    if (0 != omega)
      omega->print(std::cout);
    else
      WarningMsg << name << " not found" << endmsg;

    name = std::string("oohps");
    Particle* dummy = pTable.particle(name);
    if (0 != dummy)
      dummy->print(std::cout);
    else
      WarningMsg << name << " not found" << endmsg;

    // try clone error handling
    pTable.clone(std::string("omegaNew"), std::string("omeGa")); // syntax here is clone(new, old);
    pTable.clone(std::string("omegaNew"), std::string("omega"));
    pTable.particle(std::string("omegaNew"))->print(std::cout);

    pTable.clone(std::string("omegaNew"), std::string("omega"));

    // try to modify mass and width
    pTable.modifyMass(std::string("omegaNew"), 1.42);
    pTable.modifyWidth(std::string("omegaNew"), 0.13);
    pTable.particle(std::string("omegaNew"))->print(std::cout);

  } else {
    ErrMsg << "Error: could not parse " << pdtFile << endmsg;
    exit(1);
  }

  return 0;
}
