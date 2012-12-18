#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "PwaUtils/AbsLhNew.hh"
#include "PwaUtils/PsiToXGamStates.hh"
#include "PwaUtils/DataUtils.hh"
#include "Minuit2/MnUserParameters.h"

class AbsXdecAmp;
class pbarpReaction;
class pbarpDecAmps;

class pbarpBaseLh : public AbsLhNew {

public:
  pbarpBaseLh(boost::shared_ptr<const EvtDataBaseListNew>);

  virtual ~pbarpBaseLh();

  virtual AbsLhNew* clone_() const {
    return new  pbarpBaseLh(_evtListPtr);
  }
  virtual double calcEvtIntensity( EvtDataNew* theData, fitParamsNew& theParamVal);
  
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
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
 
  virtual complex<double> pbarpToXYAmp(Spin Minit, Spin lamX, Spin lamY,
					 EvtDataNew* theData, 
                                         double pbarpMag, double pbarpPhi);
  
  virtual void checkParamVariation(fitParamsNew& theParamVal);
  virtual void updateFitParams(fitParamsNew& theParamVal);
  virtual void cacheTheAmps(); 
private:

  // void getDefaultParams(const std::vector<std::string>& hyps, 
  // 			std::vector< boost::shared_ptr<const JPCLS> > lsAmps, 
  // 			fitParamsNew& fitVal, fitParamsNew& fitErr);
 
  void initialize();
  void fillMap(std::vector< boost::shared_ptr<const JPCLS> >& pbarpLSs, std::vector<boost::shared_ptr<pbarpDecAmps> >& decs, std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<pbarpDecAmps> >, pawian::Collection::SharedPtrLess > toFill);
  
};
