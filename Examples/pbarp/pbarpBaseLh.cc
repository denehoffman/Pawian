#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarp/pbarpBaseLh.hh"
#include "Examples/pbarp/pbarpEnv.hh"
#include "Examples/pbarp/pbarpReaction.hh"
#include "Examples/pbarp/pbarpDecAmps.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/FitParamsBaseNew.hh"

#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


pbarpBaseLh::pbarpBaseLh(boost::shared_ptr<const EvtDataBaseListNew> theEvtList) :
  AbsLhNew(theEvtList)
  ,_usePhasespace(false)
  ,_phasespaceKey("Phasespace")
{
  initialize();
}



pbarpBaseLh::~pbarpBaseLh()
{;
}

double pbarpBaseLh::calcEvtIntensity(EvtDataNew* theData, fitParamsNew& theParamVal){

  double result=0.;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > magValMap= theParamVal.Mags["pbarp"];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > phiValMap=theParamVal.Phis["pbarp"];


//   std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarDecay> >, pawian::Collection::SharedPtrLess > singletDecs=_pbarpReactionPtr->singletDecMap();
//   std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarDecay> >, pawian::Collection::SharedPtrLess >::iterator it;

//   complex<double> singletAmp(0.,0.);
//   for(it=singletDecs.begin(); it!=singletDecs.end(); ++it){
//     boost::shared_ptr<const JPCLS> theJPCLS=it->first;
//     std::vector<boost::shared_ptr<IsobarDecay> > theDecs=it->second;    
//     std::vector<boost::shared_ptr<IsobarDecay> >::iterator itDec;
//     for( itDec=theDecs.begin(); itDec!=theDecs.end(); ++itDecs){
// singletAmp+=itDecs->
//     }

//     double theMag=magValMap[*it];
//     double thePhi=phiValMap[*it];
//     complex<double> expiPhi(cos(thePhi), sin(thePhi));
//     singletAmp*=theJPCLS->preFactor*theMag*expiPhi;
//   }

  if(_usePhasespace){
    result = result + theParamVal.otherParams[_phasespaceKey];
  }
  return result;  

}

complex<double> pbarpBaseLh::pbarpToXYAmp(Spin Minit, Spin lamX, Spin lamY,
					  EvtDataNew* theData, 
					  double pbarpMag, double pbarpPhi){
  complex<double> result(0.,0.);

  return result;

}

void pbarpBaseLh::getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr){ 

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;
  
  std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
  for ( it = _jpclsStates.begin(); it!=_jpclsStates.end(); ++it){
    currentMagValMap[*it] = 1.0;
    currentPhiValMap[*it] = 0.;
    currentMagErrMap[*it] = 0.3;
    currentPhiErrMap[*it] = 0.3;
  }
  fitVal.Mags["pbarp"]=currentMagValMap;
  fitVal.Phis["pbarp"]=currentPhiValMap;
  fitErr.Mags["pbarp"]=currentMagErrMap;
  fitErr.Phis["pbarp"]=currentPhiErrMap;

  std::vector< boost::shared_ptr<pbarpDecAmps> >::iterator itDecs;
  for(itDecs=_decAmps.begin(); itDecs!=_decAmps.end(); ++itDecs){
    (*itDecs)->getDefaultParams(fitVal, fitErr);
  }
}

void pbarpBaseLh::print(std::ostream& os) const{

}


void  pbarpBaseLh::initialize(){
  boost::shared_ptr<pbarpReaction> _pbarpReactionPtr= pbarpEnv::instance()->reaction();
  _jpclsStates=_pbarpReactionPtr->jpclsStates();
  std::vector< boost::shared_ptr<IsobarDecay> > theDecs = _pbarpReactionPtr->productionDecays();

  std::vector< boost::shared_ptr<IsobarDecay> >::iterator it;
  for (it=theDecs.begin(); it!=theDecs.end(); ++it){
    boost::shared_ptr<pbarpDecAmps> currentAmp(new pbarpDecAmps(*it));
    _decAmps.push_back(currentAmp);
  }

  std::vector< boost::shared_ptr<const JPCLS> > jpclsSingleStates=_pbarpReactionPtr->jpclsSingletStates();
  fillMap(jpclsSingleStates, _decAmps, _decAmpsSinglet);

  std::vector< boost::shared_ptr<const JPCLS> > jpclsTriplet0States=_pbarpReactionPtr->jpclsTriplet0States();
  fillMap(jpclsTriplet0States, _decAmps, _decAmpsTriplet0);

  std::vector< boost::shared_ptr<const JPCLS> > jpclsTripletp1States=_pbarpReactionPtr->jpclsTripletp1States();
  fillMap(jpclsTripletp1States, _decAmps, _decAmpsTripletp1);

  std::vector< boost::shared_ptr<const JPCLS> > jpclsTripletm1States=_pbarpReactionPtr->jpclsTripletm1States();
  fillMap(jpclsTripletm1States, _decAmps, _decAmpsTripletm1);  
}

void pbarpBaseLh::fillMap(std::vector< boost::shared_ptr<const JPCLS> >& pbarpLSs, std::vector<boost::shared_ptr<pbarpDecAmps> >& decs, std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<pbarpDecAmps> >, pawian::Collection::SharedPtrLess > toFill){

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  for (itJPCLS = pbarpLSs.begin(); itJPCLS != pbarpLSs.end(); ++itJPCLS){
    std::vector<boost::shared_ptr<pbarpDecAmps> > currentAmpVec;

    std::vector<boost::shared_ptr<pbarpDecAmps> >::iterator itAmp;
    for(itAmp=decs.begin(); itAmp!=decs.end(); ++itAmp){
      if( (*(*itAmp)->jpcPtr())==(*(*itJPCLS)) )  currentAmpVec.push_back(*itAmp);
    }
    toFill[(*itJPCLS)]=currentAmpVec;   
  }

}

void pbarpBaseLh::checkParamVariation(fitParamsNew& theParamVal){
  return;
}

void pbarpBaseLh::cacheTheAmps(){
  return;
}

void pbarpBaseLh::updateFitParams(fitParamsNew& theParamVal){
}
