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

#include <vector>
#include <string>

#include <boost/config.hpp> // put this first to suppress some VC++ warnings
#include <boost/variant/recursive_variant.hpp>

class ParticleTable;

namespace decayGraph {
  struct EdgeList;
}

namespace setupGrammar {

  struct decay_tree;

  typedef
  boost::variant<
    boost::recursive_wrapper<decay_tree>
    , std::string
    >
  decay_node;

  struct decay_tree {
    std::string name;                        // mother particle name
    std::vector<decay_node> children;        // children
    std::vector<std::string> addParticle;
    std::vector<std::string> cloneParticle;
    std::vector<std::string> modParticle;
    std::vector<std::string> beamInput;
    std::vector<std::string> mcInput;
    std::vector<std::string> defineTuple;
    std::vector<std::string> fitVars;
    std::vector<std::string> initialProps;
  };

}

class SetupParser
{
public:

  SetupParser();
  ~SetupParser();

  bool parse(std::string& fileName, ParticleTable* pdtTable = 0);
  const setupGrammar::decay_tree* setup() const;
  const decayGraph::EdgeList* edgeList() const;

private:
  
  setupGrammar::decay_tree* thisDecay;
};
