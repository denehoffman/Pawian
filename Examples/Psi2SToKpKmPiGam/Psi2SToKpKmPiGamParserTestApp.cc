#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamParser.hh"
#include "Examples/Psi2SToKpKmPiGam/StreamChic1ToKpKmPiGamFitParms.hh"
#include "Examples/Psi2SToKpKmPiGam/FitParamsChic1ToKpKmPiGam.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "ErrLogger/ErrLogger.hh"

void setErrLogMode( const Psi2SToKpKmPiGamParser::enErrLogMode& erlMode ) {
  switch(erlMode) {
  case Psi2SToKpKmPiGamParser::debug :
    ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
    break;
  case Psi2SToKpKmPiGamParser::trace :
    ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
    break;
  case Psi2SToKpKmPiGamParser::routine :
    ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
    break;
  case Psi2SToKpKmPiGamParser::warning :
    ErrLogger::instance()->setLevel(log4cpp::Priority::WARN);
    break;
  case Psi2SToKpKmPiGamParser::error :
    ErrLogger::instance()->setLevel(log4cpp::Priority::ERROR);
    break;
  case Psi2SToKpKmPiGamParser::alert :
    ErrLogger::instance()->setLevel(log4cpp::Priority::ALERT);
    break;
  default: 
    ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
  }
}


int main(int __argc,char *__argv[]){

  // Parse the command line
  static Psi2SToKpKmPiGamParser theAppParams(__argc, __argv);

  // Inform the audience about the execution mode
//   emitExecutionMode(theAppParams.getAppExecMode());

   // Set the desired error logging mode
  setErrLogMode(theAppParams.getErrLogMode());

  std::string theCfgFile = theAppParams.getConfigFile();
  Info << "The path to config file is " << theCfgFile << "\n" << endmsg;    

  std::string paramStreamerPath=theAppParams.fitParamFile();

  StreamChic1ToKpKmPiGamFitParms theParamStreamer(paramStreamerPath);

  fitParams theStartparams=theParamStreamer.getFitParamVal();
  fitParams theErrorparams=theParamStreamer.getFitParamErr();    

  boost::shared_ptr<FitParamsBase> theFitParamBase=boost::shared_ptr<FitParamsBase>(new FitParamsChic1ToKpKmPiGam());

  Info << "\nThe start param values: " << "\n" << endmsg;  
  theFitParamBase->printParams(theStartparams);

  Info << "\nThe start param errors: " << "\n" << endmsg; 
  theFitParamBase->printParams(theErrorparams);


  std::ofstream theStream ( "ParameterDump.dat");
  theFitParamBase->dumpParams(theStream, theStartparams, theErrorparams);
  theStream.close();
     
  return 0;
}

