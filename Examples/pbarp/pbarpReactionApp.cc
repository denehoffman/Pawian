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
#include "PwaUtils/AbsLhNew.hh"
#include "PwaUtils/FitParamsBaseNew.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/ErrLogUtils.hh"
#include "Examples/pbarp/IsobarDecay.hh"
#include "Examples/pbarp/IsobarDecayList.hh"
#include "Examples/pbarp/pbarpEnv.hh"
#include "Examples/pbarp/pbarpReaction.hh"
#include "Examples/pbarp/pbarpDataBaseList.hh"
#include "Examples/pbarp/pbarpBaseLh.hh"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int __argc,char *__argv[]){
  clock_t start, end;
  start= clock();

  // Parse the command line
  static pbarpParser theAppParams(__argc, __argv);

  // Set the desired error logging mode
  setErrLogMode(theAppParams.getErrLogMode());
  
#ifdef _OPENMP
  const int noOfThreads=theAppParams.noOfThreads();
  omp_set_num_threads(noOfThreads);
#endif


  pbarpEnv::instance()->setup(theAppParams);

  boost::shared_ptr<pbarpReaction> thepbarpReaction=pbarpEnv::instance()->reaction();

  thepbarpReaction->print(std::cout);

  boost::shared_ptr<pbarpDataBaseList> thepbarbDataBaseListPtr(new pbarpDataBaseList()); 
  boost::shared_ptr<AbsLhNew> theLhPtr(new pbarpBaseLh(thepbarbDataBaseListPtr));

  std::string mode=theAppParams.mode();

  boost::shared_ptr<FitParamsBaseNew> theFitParamBase=boost::shared_ptr<FitParamsBaseNew>(new FitParamsBaseNew());

  if (mode=="dumpDefaultParams"){
    fitParamsNew defaultVal;
    fitParamsNew defaultErr;
    theLhPtr->getDefaultParams(defaultVal, defaultErr);
    std::string defaultparamsname = "defaultparams.dat";
    std::ofstream theStreamDefault ( defaultparamsname.c_str() );
    
    theFitParamBase->dumpParams(theStreamDefault, defaultVal, defaultErr);
    return 0;
  }

}     
 
  


