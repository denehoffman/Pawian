#include "Setup/SetupParser.hh"
#include "Setup/SetupGrammar.hh"

// Particle.hh defines ParticleData struct
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/PdtParser.hh"

#include "ErrLogger/ErrLineLog.hh"

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

bool SetupParser::parse(std::string& fileName, ParticleTable* pdtTable)
{

  std::ifstream in(fileName.c_str(), std::ios_base::in);

  if (!in) {
    ErrMsg(fatal) << "Error: Could not open input file: "
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
    ErrMsg(trace) << "\n\n"
		  << "-------------------------\n"
		  << "Parsing succeeded\n"
		  << "-------------------------\n" << endmsg;
    setupGrammar::decay_tree_printer printer;
    printer(*thisDecay);

    std::vector<std::string>::iterator cmdLine;
    ParticleData* pData;
    PdtParser pdtParser;
    for (cmdLine = thisDecay->addParticle.begin(); cmdLine != thisDecay->addParticle.end(); ++cmdLine) {
      std::cout << "add: " << *cmdLine << std::endl;
      pData = new ParticleData;
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
    ErrMsg(error) << "\n\n"
		  << "-------------------------\n"
		  << "Parsing failed\n"
		  << "stopped at: \": " << context << "...\"\n"
		  << "-------------------------\n" << endmsg;
    return false;
  }
  
  return false; // success

}

