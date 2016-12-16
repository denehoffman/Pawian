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

#include "Setup/SetupParser.hh"
#include "Setup/SetupGrammar.hh"

// Particle.hh defines ParticleData struct
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/PdtParser.hh"
#include "DecayTree/DecayTree.hh"

#include "ErrLogger/ErrLogger.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

SetupParser::SetupParser()
{
}

SetupParser::~SetupParser()
{
}


const setupGrammar::decay_tree* SetupParser::setup() const
{
  return thisDecay;
}

using namespace decayGraph;

const decayGraph::EdgeList* SetupParser::edgeList() const
{
  return &setupGrammar::edgeList;
}

bool SetupParser::parse(std::string& fileName, ParticleTable* pdtTable)
{

  std::ifstream in(fileName.c_str(), std::ios_base::in);

  if (!in) {
    Alert << "Error: Could not open input file: "
	  << fileName << endmsg;
    return 1;
  }

  std::string storage; // We will read the contents here.
  in.unsetf(std::ios::skipws); // No white space skipping!
  std::copy(std::istream_iterator<char>(in),
	    std::istream_iterator<char>(),
	    std::back_inserter(storage));
  
  typedef setupGrammar::setup_file_grammar<std::string::const_iterator> setup_file_grammar;
  setup_file_grammar setupGrammar; // Our grammar
  thisDecay = new setupGrammar::decay_tree(); // Our setup
  
  using boost::spirit::ascii::space;
  std::string::const_iterator iter = storage.begin();
  std::string::const_iterator end = storage.end();
  bool r = phrase_parse(iter, end, setupGrammar, space, *thisDecay);
  
   if (r && iter == end) {
    InfoMsg << "\n\n"
  	 << "-------------------------\n"
  	 << "Parsing succeeded\n"
  	 << "-------------------------\n" << endmsg;
    setupGrammar::decay_tree_printer printer;
    printer(*thisDecay);
    setupGrammar::edgeList.print(std::cout);

    std::vector<std::string>::iterator cmdLine;
    PdtParser pdtParser;
    for (cmdLine = thisDecay->addParticle.begin(); cmdLine != thisDecay->addParticle.end(); ++cmdLine) {
      std::cout << "add: " << *cmdLine << std::endl;
      ParticleData* pData = new ParticleData;
      if (pdtParser.parse(cmdLine->begin(), cmdLine->end(), *pData)) { // success
  	Particle* newParticle = new Particle(*pData);
  	newParticle->print(std::cout);
  	if (0 != pdtTable)
  	  pdtTable->addParticle(newParticle);
      }
    }
    if (0 != pdtTable)
      pdtTable->print(std::cout);

    return true; // true means success
   }
  
   else {
     std::string::const_iterator some = iter+30;
     std::string context(iter, (some>end)?end:some);
     ErrMsg << "\n\n"
	    << "-------------------------\n"
	    << "Parsing failed\n"
	    << "stopped at: \": " << context << "...\"\n"
	    << "-------------------------\n" << endmsg;
     return false;
   }
   
   return false; // success

}

