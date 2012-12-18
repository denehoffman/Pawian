#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarp/pbarpDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/DataUtils.hh"
#include "Examples/pbarp/IsobarDecay.hh"

#ifdef _OPENMP
#include <omp.h>
#endif

pbarpDecAmps::pbarpDecAmps(boost::shared_ptr<IsobarDecay> theDec) :
  AbsXdecAmp(theDec->name())
  ,_decay(theDec)
  ,_JPCPtr(theDec->motherJPC())
  ,_JPCLSs(theDec->JPCLSAmps())
  ,_key("_"+theDec->fitParSuffix())
  ,_daughter1IsStable(theDec->isDaughter1Stable())
  ,_daughter2IsStable(theDec->isDaughter2Stable())
{
  initialize();
}

pbarpDecAmps::~pbarpDecAmps()
{
}

complex<double> pbarpDecAmps::XdecAmp(Spin lamX, EvtDataNew* theData){

  complex<double> result(0.,0.);

  return result;
}

void  pbarpDecAmps::getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr){

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
  for(itLS=_JPCLSs.begin(); itLS!=_JPCLSs.end(); ++itLS){
    currentMagValMap[*itLS]=0.5;
    currentPhiValMap[*itLS]=0.;
    currentMagErrMap[*itLS]=0.5;
    currentPhiErrMap[*itLS]=0.3;
  }

  fitVal.Mags[_key]=currentMagValMap;
  fitVal.Phis[_key]=currentPhiValMap;
  fitErr.Mags[_key]=currentMagErrMap;
  fitErr.Phis[_key]=currentPhiErrMap;

  if(!_daughter1IsStable) _decAmpDaughter1->getDefaultParams(fitVal, fitErr);
  if(!_daughter2IsStable) _decAmpDaughter2->getDefaultParams(fitVal, fitErr);  
}

void pbarpDecAmps::print(std::ostream& os) const{
  return; //dummy
}

void pbarpDecAmps::initialize(){
  if(!_daughter1IsStable){
    boost::shared_ptr<IsobarDecay> decDaughter1=_decay->decDaughter1();
    _decAmpDaughter1=boost::shared_ptr<pbarpDecAmps>(new pbarpDecAmps(decDaughter1));
  }

  if(!_daughter2IsStable){
    boost::shared_ptr<IsobarDecay> decDaughter2=_decay->decDaughter1();
    _decAmpDaughter2=boost::shared_ptr<pbarpDecAmps>(new pbarpDecAmps(decDaughter2));
  }

}

void pbarpDecAmps::checkRecalculation(fitParamsNew& theParamVal){
}
 

void  pbarpDecAmps::updateFitParams(fitParamsNew& theParamVal){
}
