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

#include "Particle/PdtParser.hh"

#include "Particle/ParticleGrammar.hh"

PdtParser::PdtParser()
{
}

PdtParser::~PdtParser()
{
}

bool PdtParser::parse(std::string& fileName, ParticleTable& table)
{
  ParticleData* pData;

  std::ifstream file;
  file.open(fileName.c_str());

  if (file.is_open()) {
    std::string str;
    while (!file.eof()) {
      std::getline(file, str);
      if (!file.eof()) {
	pData = new ParticleData;
	if (particleGrammar::parse_particle(str.begin(), str.end(), *pData)) { // success
	  Particle* newParticle = new Particle(*pData);
	  table.addParticle(newParticle);
	}
      }
    }
    file.close();
    return true; // success
  }
  else
    return false;
}

bool PdtParser::parse(std::string::const_iterator begin,
		      std::string::const_iterator end,
		      ParticleData& pData)
{
  return particleGrammar::parse_particle(begin, end, pData);
}
