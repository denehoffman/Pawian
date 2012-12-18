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

  std::map <boost::shared_ptr<const JPCLS>, std::vector< boost::shared_ptr<pbarpDecAmps> >, pawian::Collection::SharedPtrLess >::iterator it;

  //  std::cout <<"_decAmpsSinglet.size():\t" << _decAmpsSinglet.size() << std::endl;

  complex<double> singletAmp(0.,0.);
  for(it=_decAmpsSinglet.begin(); it!=_decAmpsSinglet.end(); ++it){
    complex<double> tmpAmp(0.,0.);
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;
    std::vector<boost::shared_ptr<pbarpDecAmps> > decAmps=it->second;    
     std::vector<boost::shared_ptr<pbarpDecAmps> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       complex<double> currentDecAmp=(*itDec)->XdecAmp(0, theData);
       // std::cout << "decAmp: " << currentDecAmp << std::endl;
       tmpAmp+=currentDecAmp;
     }

     double theMag=_currentParamMags[theJPCLS];
     double thePhi=_currentParamPhis[theJPCLS];
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=theJPCLS->preFactor*theMag*expi;
     singletAmp+=tmpAmp;
  }
  // std::cout << "singletAmp: " << singletAmp << std::endl;

  complex<double> triplet0Amp(0.,0.);
  for(it=_decAmpsTriplet0.begin(); it!=_decAmpsTriplet0.end(); ++it){
    complex<double> tmpAmp(0.,0.);
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;
    std::vector<boost::shared_ptr<pbarpDecAmps> > decAmps=it->second;    
     std::vector<boost::shared_ptr<pbarpDecAmps> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       tmpAmp+=(*itDec)->XdecAmp(0, theData);
     }

     double theMag=_currentParamMags[theJPCLS];
     double thePhi=_currentParamPhis[theJPCLS];
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=theJPCLS->preFactor*theMag*expi;
     triplet0Amp+=tmpAmp;
  }

  complex<double> tripletp1Amp(0.,0.);
  for(it=_decAmpsTripletp1.begin(); it!=_decAmpsTripletp1.end(); ++it){
    complex<double> tmpAmp(0.,0.);
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;
    std::vector<boost::shared_ptr<pbarpDecAmps> > decAmps=it->second;    
     std::vector<boost::shared_ptr<pbarpDecAmps> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       tmpAmp+=(*itDec)->XdecAmp(1, theData);
     }

     double theMag=_currentParamMags[theJPCLS];
     double thePhi=_currentParamPhis[theJPCLS];
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=theJPCLS->preFactor*theMag*expi;
     tripletp1Amp+=tmpAmp;
  }  

  complex<double> tripletm1Amp(0.,0.);
  for(it=_decAmpsTripletm1.begin(); it!=_decAmpsTripletm1.end(); ++it){
    complex<double> tmpAmp(0.,0.);
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;
    std::vector<boost::shared_ptr<pbarpDecAmps> > decAmps=it->second;    
     std::vector<boost::shared_ptr<pbarpDecAmps> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       tmpAmp+=(*itDec)->XdecAmp(-1, theData);
     }

     double theMag=_currentParamMags[theJPCLS];
     double thePhi=_currentParamPhis[theJPCLS];
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=theJPCLS->preFactor*theMag*expi;
     tripletm1Amp+=tmpAmp;
  } 


  result = norm(singletAmp)+ norm(triplet0Amp)+ norm(tripletp1Amp)+ norm(tripletm1Amp);

  if(_usePhasespace){
    result = result + theParamVal.otherParams[_phasespaceKey];
  }
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

  std::vector< boost::shared_ptr<const JPCLS> > jpclsSingletStates=_pbarpReactionPtr->jpclsSingletStates();
  fillMap(jpclsSingletStates, _decAmps, _decAmpsSinglet);

  std::vector< boost::shared_ptr<const JPCLS> > jpclsTriplet0States=_pbarpReactionPtr->jpclsTriplet0States();
  fillMap(jpclsTriplet0States, _decAmps, _decAmpsTriplet0);

  std::vector< boost::shared_ptr<const JPCLS> > jpclsTripletp1States=_pbarpReactionPtr->jpclsTripletp1States();
  fillMap(jpclsTripletp1States, _decAmps, _decAmpsTripletp1);

  std::vector< boost::shared_ptr<const JPCLS> > jpclsTripletm1States=_pbarpReactionPtr->jpclsTripletm1States();
  fillMap(jpclsTripletm1States, _decAmps, _decAmpsTripletm1);  
}

void pbarpBaseLh::fillMap(std::vector< boost::shared_ptr<const JPCLS> >& pbarpLSs, std::vector<boost::shared_ptr<pbarpDecAmps> >& decs, std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<pbarpDecAmps> >, pawian::Collection::SharedPtrLess >& toFill){


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
  std::vector< boost::shared_ptr<pbarpDecAmps> >::iterator it;
  for (it=_decAmps.begin(); it!=_decAmps.end(); ++it){
    (*it)->checkRecalculation(theParamVal);
  }
}

void pbarpBaseLh::cacheTheAmps(){
  std::vector< boost::shared_ptr<pbarpDecAmps> >::iterator it;
  for (it=_decAmps.begin(); it!=_decAmps.end(); ++it){
    (*it)->cacheAmplitudes();
  }
}

void pbarpBaseLh::updateFitParams(fitParamsNew& theParamVal){

   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > magMap=theParamVal.Mags["pbarp"];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > phiMap=theParamVal.Phis["pbarp"];

   std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
   for (it=_jpclsStates.begin(); it!=_jpclsStates.end(); ++it){
     double theMag=magMap[*it];
     double thePhi=phiMap[*it];
     _currentParamMags[*it]=theMag;
     _currentParamPhis[*it]=thePhi;
   }

  std::vector< boost::shared_ptr<pbarpDecAmps> >::iterator itDecs;
  for(itDecs=_decAmps.begin(); itDecs!=_decAmps.end(); ++itDecs){
    (*itDecs)->updateFitParams(theParamVal);
  }
}
