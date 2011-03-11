/**
 * @file GArgumentParser.hpp
 */

/* Copyright (C) Dr. Ruediger Berlich and Karlsruhe Institute of Technology
 * (University of the State of Baden-Wuerttemberg and National Laboratory
 * of the Helmholtz Association)
 *
 * Contact: info [at] gemfony (dot) com
 *
 * This file is part of the Geneva library, Gemfony scientific's optimization
 * library.
 *
 * Geneva is free software: you can redistribute it and/or modify
 * it under the terms of version 3 of the GNU Affero General Public License
 * as published by the Free Software Foundation.
 *
 * Geneva is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with the Geneva library.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For further information on Gemfony scientific and Geneva, visit
 * http://www.gemfony.com .
 */

// Standard headers go here

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Boost headers go here

#include <boost/version.hpp>

#if BOOST_VERSION < 103600
#error "Error: Boost should at least have version 1.36 !"
#endif /* BOOST_VERSION */

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#ifndef GMARGUMENTPARSER_HPP_
#define GMARGUMENTPARSER_HPP_

// For Microsoft-compatible compilers
#if defined(_MSC_VER)  &&  (_MSC_VER >= 1020)
#pragma once
#endif

// GenEvA headers go here
#include <common/GCommonEnums.hpp>
#include <common/GSerializationHelperFunctionsT.hpp>
#include <geneva/GOptimizationEnums.hpp>

namespace po = boost::program_options;
using namespace Gem::Geneva;

class ApplicationParameter
{
  public:
  typedef enum tagExecMode {GenEvA=0, Minuit=1, GenToMinuit=2, SpinDensity=3, QAmode=4} enExecMode;
    typedef enum tagerrLogMode {debug,trace,routine,warning,error,alert} enErrLogMode;

  public:
    ApplicationParameter(int argc,char **argv):configFile("./GOmegaPiProject.cfg"),
					       parallelizationMode(1),
                                               ip("localhost"),
					       port(10000),
					       serMode(Gem::Common::SERIALIZATIONMODE_TEXT),
					       enAppExecMode(GenEvA),
					       nProducerThreads(10),
					       nEvaluationThreads(4),
					       populationSize(100),
					       nParents(5),
					       maxIterations(2000),
					       maxMinutes(10),
					       reportIteration(1),
					       smode(MUPLUSNU),
					       arraySize(1000),
					       iM(0),
                                               iM_(0),
                                               bCalcAllSpindensity(false),
                                               bRemoveEventsFromEnd(true),
                                               nNumEventsRemove(0),
                                               nNumEventsRed(0),
					       processingCycles(1),
					       returnRegardless(true),
					       waitFactor(5),
					       jMax(3),
					       pbarMom(600),
					       errLogMode(debug),
					       theSourcePath("./")
  {
    if (!parseCommandLine(argc, argv)) throw false;
  }
  

    inline const std::string& getConfigFile() { return configFile;}
    inline const boost::uint16_t& getParallelizationMode() {return parallelizationMode;}
    inline const bool& getServerMode() { return serverMode; }
    inline const std::string& getIp() { return ip; }
    inline const unsigned short& getPort() { return port; }
    inline const boost::uint16_t& getNProducerThreads() { return nProducerThreads; }
    inline const boost::uint16_t& getNEvaluationThreads() { return nEvaluationThreads; }
    inline const std::size_t& getPopulationSize() { return populationSize; }
    inline const std::size_t& getNParents() { return nParents; }
    inline const boost::uint32_t& getMaxIterations() { return maxIterations; }
    inline const long& getMaxMinutes() { return maxMinutes; }
    inline const boost::uint32_t& getReportIteration() { return reportIteration; }
    inline const recoScheme& getRScheme() { return rScheme; }
    inline const std::size_t& getArraySize() { return arraySize; }
    inline const sortingMode& getSmode() { return smode; }
    inline const boost::uint32_t& getProcessingCycles() { return processingCycles; }
    inline const bool& getReturnRegardless() { return returnRegardless; }
    inline const boost::uint32_t& getWaitFactor() { return waitFactor; }
    inline const unsigned& getJMax() { return jMax; }
    inline const unsigned& getPbarMom() { return pbarMom; }
    inline const enErrLogMode& getErrLogMode() { return errLogMode; }
    inline const Gem::Common::serializationMode& getSerMode() { return serMode; }
    inline const enExecMode& getAppExecMode() { return enAppExecMode; }
    inline const std::string& getSourcePath() { return theSourcePath; }
    inline const std::string& getPathStartParamFile() { return strPathStartParamFile; }
    inline const std::string& getName() { return strName; }
    inline const int getM() { return iM; }
    inline const int getM_() { return iM_; }
    inline const bool getCalcAllSpindensity() { return bCalcAllSpindensity; }
    inline const bool getRemoveEventsFromEnd() { return bRemoveEventsFromEnd; }
    inline const int getNumEventsRemove() { return nNumEventsRemove; }
    inline const int getNumEventsRed() { return nNumEventsRed; }
    
  protected:
    bool parseCommandLine(int argc,char **argv);

  protected:
  //Actual parsed settings
  std::string strName;
  std::string strPathStartParamFile;
  std::string configFile;		  
  boost::uint16_t parallelizationMode;
  bool serverMode;
  std::string ip;
  unsigned short port;
  Gem::Common::serializationMode serMode;
  enExecMode enAppExecMode;
  boost::uint16_t nProducerThreads;
  boost::uint16_t nEvaluationThreads;
  std::size_t populationSize;
  std::size_t nParents;
  boost::uint32_t maxIterations;
  long maxMinutes;
  boost::uint32_t reportIteration;
  recoScheme rScheme;
  sortingMode smode;
  std::size_t arraySize;
  int iM;
  int iM_;
  
  bool bCalcAllSpindensity;
  bool bRemoveEventsFromEnd;
  int nNumEventsRemove;
  int nNumEventsRed;
  boost::uint32_t processingCycles;
  bool returnRegardless;
  boost::uint32_t waitFactor;
  unsigned jMax;
  unsigned pbarMom;
  enErrLogMode errLogMode;
  std::string theSourcePath;
};
#endif /* GARGUMENTPARSER_HPP_ */
