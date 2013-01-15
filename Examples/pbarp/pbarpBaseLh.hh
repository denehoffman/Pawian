#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/PsiToXGamStates.hh"
#include "PwaUtils/DataUtils.hh"
#include "Minuit2/MnUserParameters.h"

class AbsXdecAmp;
class pbarpReaction;
class pbarpDecAmps;

class pbarpBaseLh : public AbsLh {

public:
  pbarpBaseLh(boost::shared_ptr<const EvtDataBaseList>);

  virtual ~pbarpBaseLh();

  virtual AbsLh* clone_() const {
    return new  pbarpBaseLh(_evtListPtr);
  }
  virtual double calcEvtIntensity( EvtData* theData, fitParams& theParamVal);
  
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual void print(std::ostream& os) const;

protected:
   bool _usePhasespace;  
  const std::string _phasespaceKey;
  boost::shared_ptr<pbarpReaction> _pbarpReactionPtr;
  std::vector< boost::shared_ptr<const JPCLS> > _jpclsStates;
  std::vector< boost::shared_ptr<pbarpDecAmps> > _decAmps;
  std::map <boost::shared_ptr<const JPCLS>, std::vector< boost::shared_ptr<pbarpDecAmps> >, pawian::Collection::SharedPtrLess > _decAmpsSinglet;
  std::map <boost::shared_ptr<const JPCLS>, std::vector< boost::shared_ptr<pbarpDecAmps> >, pawian::Collection::SharedPtrLess > _decAmpsTriplet0;
  std::map <boost::shared_ptr<const JPCLS>, std::vector< boost::shared_ptr<pbarpDecAmps> >, pawian::Collection::SharedPtrLess > _decAmpsTripletp1;
  std::map <boost::shared_ptr<const JPCLS>, std::vector< boost::shared_ptr<pbarpDecAmps> >, pawian::Collection::SharedPtrLess > _decAmpsTripletm1;
 
 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > _currentParamMags;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > _currentParamPhis;
  
  virtual void checkParamVariation(fitParams& theParamVal);
  virtual void updateFitParams(fitParams& theParamVal);
  virtual void cacheTheAmps();
 
private:

  // void getDefaultParams(const std::vector<std::string>& hyps, 
  // 			std::vector< boost::shared_ptr<const JPCLS> > lsAmps, 
  // 			fitParams& fitVal, fitParams& fitErr);
 
  void initialize();
  void fillMap(std::vector< boost::shared_ptr<const JPCLS> >& pbarpLSs, std::vector<boost::shared_ptr<pbarpDecAmps> >& decs, std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<pbarpDecAmps> >, pawian::Collection::SharedPtrLess >& toFill);
  
};
