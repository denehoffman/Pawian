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

// pbarpParser class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include "PwaUtils/ParserBase.hh"
// Boost headers go here



class pbarpParser : public ParserBase 
{

  public:

  pbarpParser(int argc,char **argv);
  virtual ~pbarpParser(){;}

  const unsigned int getLMax() const { return _lMax; }
  const float getpbarMomentum() const { return _pbarMomentum; } 
  const std::vector<std::string>& spinDensityNames() const { return _spinDensity;}
  const std::vector<std::string>& lmaxParticle() const{ return _lmaxParticle; }
  const std::vector<std::string>& dropPbarpLForParticle() const{ return _dropPbarpLForParticle; }

protected:
  virtual bool parseCommandLine(int argc,char **argv); 
  unsigned int _lMax;
  float _pbarMomentum;
  std::vector<std::string> _spinDensity;
  std::vector<std::string> _lmaxParticle;
  std::vector<std::string> _dropPbarpLForParticle;
};



