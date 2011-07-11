/**
 * @file GArgumentParserLS.hpp
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

#ifndef GMARGUMENTPARSERLS_HPP_
#define GMARGUMENTPARSERLS_HPP_

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

class ApplicationParameterLS
{
  public:
  typedef enum tagExecMode {GenEvA=0, Minuit=1, GenToMinuit=2, SpinDensity=3, QAmode=4, StatesTest=5, HistTest=6, ParamStreamTest=7} enExecMode;
    typedef enum tagerrLogMode {debug,trace,routine,warning,error,alert} enErrLogMode;

  public:
    ApplicationParameterLS(int argc,char **argv)
      : configFile("./GOmegaPiProject.cfg")
      ,	parallelizationMode(1)
      , ip("localhost")
      , port(10000)
      , serMode(Gem::Common::SERIALIZATIONMODE_TEXT)
      , enAppExecMode(GenEvA)
      , nProducerThreads(10)
      , nEvaluationThreads(4)
      , populationSize(100)
      , nParents(5)
      , maxIterations(2000)
      , maxMinutes(10)
      , reportIteration(1)
      , smode(MUPLUSNU)
      , arraySize(1000)
      , iM(0)
      , iM_(0)
      , bCalcAllSpindensity(false)
      , bRemoveEventsFromEnd(true)
      , nNumEventsRemove(0)
      , nNumEventsRed(0)
      , processingCycles(1)
      , returnRegardless(true)
      , waitFactor(5)
      , lMax(3)
      , pbarMom(600)
      , errLogMode(debug)
      , theSourcePath("./")
      , theLhMode("OmegaPiLhGammaBw")
  {
    if (!parseCommandLine(argc, argv)) throw false;
  }
  
  const std::string& getConfigFile() const { return configFile;}
  const boost::uint16_t& getParallelizationMode() const {return parallelizationMode;}
  const bool& getServerMode() const { return serverMode; }
  const std::string& getIp() const { return ip; }
  const unsigned short& getPort() const { return port; }
  const boost::uint16_t& getNProducerThreads() const { return nProducerThreads; }
  const boost::uint16_t& getNEvaluationThreads() const { return nEvaluationThreads; }
  const std::size_t& getPopulationSize() const { return populationSize; }
  const std::size_t& getNParents() const { return nParents; }
  const boost::uint32_t& getMaxIterations() const { return maxIterations; }
  const long& getMaxMinutes() const { return maxMinutes; }
  const boost::uint32_t& getReportIteration() const { return reportIteration; }
  const recoScheme& getRScheme() const { return rScheme; }
  const std::size_t& getArraySize() const { return arraySize; }
  const sortingMode& getSmode() const { return smode; }
  const boost::uint32_t& getProcessingCycles() const { return processingCycles; }
  const bool& getReturnRegardless() const { return returnRegardless; }
  const boost::uint32_t& getWaitFactor() const { return waitFactor; }
  const unsigned& getLMax() const { return lMax; }
  const unsigned& getPbarMom() const { return pbarMom; }
  const enErrLogMode& getErrLogMode() const { return errLogMode; }
  const Gem::Common::serializationMode& getSerMode() const { return serMode; }
  const enExecMode& getAppExecMode() const { return enAppExecMode; }
  const std::string& getSourcePath() const { return theSourcePath; }
  const std::string& getLhMode() const { return theLhMode; }
  const std::string& getPathStartParamFile() const { return strPathStartParamFile; }
  const std::string& getName() const { return strName; }
  const int getM() const { return iM; }
  const int getM_() const { return iM_; }
  const bool getCalcAllSpindensity() const { return bCalcAllSpindensity; }
  const bool getRemoveEventsFromEnd() const { return bRemoveEventsFromEnd; }
  const int getNumEventsRemove() const { return nNumEventsRemove; }
  const int getNumEventsRed() const { return nNumEventsRed; }
    
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
  unsigned lMax;
  unsigned pbarMom;
  enErrLogMode errLogMode;
  std::string theSourcePath;
  std::string theLhMode;
};

#endif /* GARGUMENTPARSER_HPP_ */
