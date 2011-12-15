#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKFitParams.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKStates.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiGamKsKlKKFitParams::JpsiGamKsKlKKFitParams() 
{
  JpsiGamKsKlKKStates theStates;
  theStates.print(std::cout);
  
  _jpclsMap[paramEnumGamKsKlKK::etacGamma]=theStates.PsiToEtaGammaStates();
  _jpclsMap[paramEnumGamKsKlKK::eta2225Gamma]=theStates.PsiToEtaGammaStates();
  _jpclsMap[paramEnumGamKsKlKK::f02020Gamma]=theStates.PsiToF0GammaStates();
  _jpclsMap[paramEnumGamKsKlKK::f22010Gamma]=theStates.PsiToF2GammaStates();
  _jpclsMap[paramEnumGamKsKlKK::f22300Gamma]=theStates.PsiToF2GammaStates();
  _jpclsMap[paramEnumGamKsKlKK::f22340Gamma]=theStates.PsiToF2GammaStates(); 
  _jpclsMap[paramEnumGamKsKlKK::f02020ToPhiPhi]=theStates.F0ToPhiPhiStates();


}

JpsiGamKsKlKKFitParams::~JpsiGamKsKlKKFitParams()
{
}

std::vector< boost::shared_ptr<const JPCLS> > JpsiGamKsKlKKFitParams::jpclsVec(unsigned int index){
  unsigned int nOfAmps=paramEnumGamKsKlKK::nAmps;
  if ( index > paramEnumGamKsKlKK::nAmps ) Info << "jpclsVec: index " << index 
						 << "  number of amplitude parametes paramEnumGamKsKlKK::numAmps = " 
 						    << nOfAmps ;  // << endmsg;

  return _jpclsMap[index]; 

}

std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& JpsiGamKsKlKKFitParams::ampMap(paramGamKsKlKK& params, unsigned int index){
  unsigned int nOfAmps=paramEnumGamKsKlKK::nAmps;
  if ( index > paramEnumGamKsKlKK::nAmps ) Info << "ampMap: index " << index 
						 << "  number of amplitude parametes paramEnumGamKsKlKK::numAmps = " 
 						    << nOfAmps ;  // << endmsg;

  if (index==paramEnumGamKsKlKK::etacGamma) return params.PsiToEtacGamma;
  else if (index==paramEnumGamKsKlKK::eta2225Gamma) return params.PsiToEta2225Gamma;
  else if (index==paramEnumGamKsKlKK::f02020Gamma) return params.PsiToF02020Gamma;
  else if (index==paramEnumGamKsKlKK::f22010Gamma) return params.PsiToF22010Gamma;
  else if (index==paramEnumGamKsKlKK::f22300Gamma) return params.PsiToF22300Gamma;
  else if (index==paramEnumGamKsKlKK::f22340Gamma) return params.PsiToF22340Gamma;
  else if (index==paramEnumGamKsKlKK::f02020ToPhiPhi) return params.F02020ToPhiPhi;
  else Info << "index " << index << " not supported !!!" ;  // << endmsg;

  return params.PsiToEtacGamma; //makes the compiler happy
}

pair<double, double>& JpsiGamKsKlKKFitParams::massPair(paramGamKsKlKK& params, unsigned int index){
  unsigned int nOfAmps=paramEnumGamKsKlKK::nAmps;
  unsigned int nOfMasses=paramEnumGamKsKlKK::nMasses;
  if ( index<paramEnumGamKsKlKK::nAmps || index > paramEnumGamKsKlKK::nMasses ) Info << "index " << index 
										      << "  must be between  " << nOfAmps
										      << " and " << nOfMasses ;  // << endmsg;
  if (index==paramEnumGamKsKlKK::etac) return params.BwEtac;
  else if (index==paramEnumGamKsKlKK::eta2225) return params.BwEta2225;
  else if (index==paramEnumGamKsKlKK::f02020) return params.BwF02020;
  else if (index==paramEnumGamKsKlKK::f22010) return params.BwF22010;
  else if (index==paramEnumGamKsKlKK::f22300) return params.BwF22300;
  else if (index==paramEnumGamKsKlKK::f22340) return params.BwF22340;
  else Info << "index " << index << " not supported !!!" ;  // << endmsg;

  return params.BwEtac; //makes the compiler happy
}

int JpsiGamKsKlKKFitParams::setFitParamValDec(paramGamKsKlKK& theParamVal, const std::vector<double>& par, int counter, unsigned int index){

  int resultCount=counter;
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > currentStates;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >* currentMap=0;

  currentStates = _jpclsMap[index];
  currentMap = &ampMap(theParamVal, index);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;

  for ( itJPCLS=currentStates.begin(); itJPCLS!=currentStates.end(); ++itJPCLS){
    double mag=par[resultCount];
    resultCount++;
    
    double phi=par[resultCount];
    resultCount++;
    

    std::pair <double,double> tmpParameter=make_pair(mag,phi);
    (*currentMap)[(*itJPCLS)]=tmpParameter;
  }

  return resultCount;
}


int JpsiGamKsKlKKFitParams::setFitParamValMass(paramGamKsKlKK& theParamVal, const std::vector<double>& par, int counter, unsigned int index){

  int resultCount=counter;

  pair<double, double>* currentPair= & massPair(theParamVal, index);


  double currentMass=par[resultCount];
  resultCount++;
  double currentWidth=par[resultCount];
  resultCount++;

  std::pair <double,double> thePair=make_pair(currentMass, currentWidth);

  (*currentPair)=thePair;

  return resultCount;
}



void JpsiGamKsKlKKFitParams::setMnUsrParamsDec(MnUserParameters& upar, paramGamKsKlKK& startVal,  paramGamKsKlKK& errVal, unsigned int index){
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > startParams;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > errParams;

  startParams= ampMap(startVal, index);
  errParams= ampMap(errVal, index);


  int counter=0;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for (it=startParams.begin(); it!=startParams.end(); ++it){
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;   
    std::pair <double,double> startPair=it->second;
    std::pair <double,double> errPair=errParams.find(theJPCLS)->second;

    //now fill the fitParameterMap
    std::string magStr=theJPCLS->name()+paramEnumGamKsKlKK::name(index)+"Mag";
    std::string phiStr=theJPCLS->name()+paramEnumGamKsKlKK::name(index)+"phi";

    double magVal=startPair.first;
    double phiVal=startPair.second;
    
    double magErr=errPair.first;
    double phiErr=errPair.second;
    
    double magMin=magVal-3.*magErr;
    if (magMin<0.) magMin=0.;
    
    upar.Add(magStr, magVal, magErr, magMin, magVal+3.*magErr);
    upar.Add(phiStr, phiVal, phiErr, -3.*M_PI, 3.*M_PI);

    counter++;
  }

}


void JpsiGamKsKlKKFitParams::setMnUsrParamsMass(MnUserParameters& upar, paramGamKsKlKK& startVal,  paramGamKsKlKK& errVal, unsigned int index){
  pair<double, double> startParams;
  pair<double, double> errParams;
  startParams=massPair(startVal, index);
  errParams=massPair(errVal, index);

  //now fill the fitParameterMap
  std::string massStr=paramEnumGamKsKlKK::name(index)+"Mag";
  std::string widthStr=paramEnumGamKsKlKK::name(index)+"Width";
  
  double massVal=startParams.first;
  double widthVal=startParams.second;
  
  double massErr=errParams.first;
  double widthErr=errParams.second;

  double massMin=massVal-3.*widthVal;
  if (massMin<0.) massMin=0.;

  double massMax=massVal+3.*widthVal;

  double widthMin=0.; 
  double widthMax=2.*widthVal;
 
  upar.Add(massStr, massVal, massErr, massMin, massMax);
  upar.Add(widthStr, widthVal, widthErr, widthMin, widthMax);
}



