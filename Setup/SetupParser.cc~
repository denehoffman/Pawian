#include "Setup/SetupParser.hh"

// Particle.hh defines ParticleData struct
#include "Particle/Particle.hh"

#include "Particle/ParticleTable.hh"

// SetupData.hh has all the boost includes
#include "SetupData.hh"

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



bool SetupParser::parse(std::string& fileName)
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
  
  typedef setupData::setup_file_grammar<std::string::const_iterator> setup_file_grammar;
  setup_file_grammar setupGrammar; // Our grammar
  setupData::decay_tree ast; // Our setup
  
  using boost::spirit::ascii::space;
  std::string::const_iterator iter = storage.begin();
  std::string::const_iterator end = storage.end();
  bool r = phrase_parse(iter, end, setupGrammar, space, ast);
  
  if (r && iter == end) {
    ErrMsg(trace) << "\n\n"
		  << "-------------------------\n"
		  << "Parsing succeeded\n"
		  << "-------------------------\n" << endmsg;
    setupData::decay_tree_printer printer;
    printer(ast);
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

