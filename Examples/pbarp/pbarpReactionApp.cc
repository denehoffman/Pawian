#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <map>
#include <iterator>

#include "Examples/pbarp/pbarpParser.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Particle/PdtParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/pbarpStatesLS.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Examples/pbarp/IsobarDecay.hh"
#include "Examples/pbarp/IsobarDecayList.hh"
#include "Examples/pbarp/pbarpEnv.hh"
#include "Examples/pbarp/pbarpReaction.hh"

int main(int __argc,char *__argv[]){
  clock_t start, end;
  start= clock();
  
  // Parse the command line
  static pbarpParser theAppParams(__argc, __argv);

  pbarpEnv::instance()->setup(theAppParams);

  boost::shared_ptr<pbarpReaction> thepbarpReaction=pbarpEnv::instance()->reaction();

  thepbarpReaction->print(std::cout);  

}     
 
  


