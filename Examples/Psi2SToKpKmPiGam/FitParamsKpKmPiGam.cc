#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/FitParamsKpKmPiGam.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamStates.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

FitParamsKpKmPiGam::FitParamsKpKmPiGam() 
{
  Psi2SToKpKmPiGamStates theStates;
  theStates.print(std::cout);

  _jpclsMap[paramEnumKpKmPiGam::ChiGam]=theStates.PsiToChiGamStates();
  _jpclsMap[paramEnumKpKmPiGam::K890K]=theStates.ChiToKst1KStates();
  _jpclsMap[paramEnumKpKmPiGam::K_0_1400K]=theStates.ChiToKst0KStates();
  _jpclsMap[paramEnumKpKmPiGam::K_1_1400K]=theStates.ChiToKst1KStates();
  _jpclsMap[paramEnumKpKmPiGam::K_2_1400K]=theStates.ChiToKst2KStates();
  _jpclsMap[paramEnumKpKmPiGam::a980Pi]=theStates.ChiToa0PiStates();
  _jpclsMap[paramEnumKpKmPiGam::KappaK]=theStates.ChiToKst0KStates();
  _jpclsMap[paramEnumKpKmPiGam::KKPi]=theStates.ChiToa0PiStates();
}

FitParamsKpKmPiGam::~FitParamsKpKmPiGam()
{
}

std::vector< boost::shared_ptr<const JPCLS> > FitParamsKpKmPiGam::jpclsVec(unsigned int index){
  unsigned int nOfAmps=paramEnumKpKmPiGam::nAmps;
  if ( index > paramEnumKpKmPiGam::nAmps ) Alert << "index " << index 
						 << "  number of amplitude parametes paramEnumKpKmPiGam::numAmps = " 
 						    << nOfAmps << endmsg;

  return _jpclsMap[index]; 

}

std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& FitParamsKpKmPiGam::ampMap(paramKpKmPiGam& params, unsigned int index){
  unsigned int nOfAmps=paramEnumKpKmPiGam::nAmps;
  if ( index > paramEnumKpKmPiGam::nAmps ) Alert << "index " << index 
						 << "  number of amplitude parametes paramEnumKpKmPiGam::numAmps = " 
 						    << nOfAmps << endmsg;

  if (index==paramEnumKpKmPiGam::ChiGam) return params.PsiToChiGam;
  else if (index==paramEnumKpKmPiGam::K890K) return params.ChiToK890K;
  else if (index==paramEnumKpKmPiGam::K_0_1400K) return params.ChiToK1400_0_K;
  else if (index==paramEnumKpKmPiGam::K_1_1400K) return params.ChiToK1400_1_K;
  else if (index==paramEnumKpKmPiGam::K_2_1400K) return params.ChiToK1400_2_K;
  else if (index==paramEnumKpKmPiGam::a980Pi) return params.ChiToa0Pi;
  else if (index==paramEnumKpKmPiGam::KappaK) return params.ChiToKappa_K;
  else if (index==paramEnumKpKmPiGam::KKPi) return params.ChiToKKPi;
  else Alert << "index " << index << " not supported !!!" << endmsg;

  return params.PsiToChiGam; //makes the compiler happy
}

pair<double, double>& FitParamsKpKmPiGam::massPair(paramKpKmPiGam& params, unsigned int index){
  unsigned int nOfAmps=paramEnumKpKmPiGam::nAmps;
  unsigned int nOfMasses=paramEnumKpKmPiGam::nMasses;
  if ( index<paramEnumKpKmPiGam::nAmps || index > paramEnumKpKmPiGam::nMasses ) Alert << "index " << index 
										      << "  must be between  " << nOfAmps
										      << " and " << nOfMasses << endmsg;
  if (index==paramEnumKpKmPiGam::K890) return params.BwK890;
  else if (index==paramEnumKpKmPiGam::K_0_1400) return params.BwK1400_0;
  else if (index==paramEnumKpKmPiGam::K_1_1400) return params.BwK1400_1;
  else if (index==paramEnumKpKmPiGam::K_2_1400) return params.BwK1400_2;
  else if (index==paramEnumKpKmPiGam::Kappa) return params.BwKappa;
  else Alert << "index " << index << " not supported !!!" << endmsg;

  return params.BwK1400_0; //makes the compiler happy
}

int FitParamsKpKmPiGam::setFitParamValDec(paramKpKmPiGam& theParamVal, const std::vector<double>& par, int counter, unsigned int index){

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


int FitParamsKpKmPiGam::setFitParamValMass(paramKpKmPiGam& theParamVal, const std::vector<double>& par, int counter, unsigned int index){

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


int FitParamsKpKmPiGam::setFitParamFlattea980Mass(paramKpKmPiGam& theParamVal, const std::vector<double>& par, int counter, std::string key){

  int resultCount=counter;

  if (key==paramEnumKpKmPiGam::name(paramEnumKpKmPiGam::a980)){ 
    double* resultFlatMa980= &theParamVal.FlatMa980;
    double* resultFlatgKK= &theParamVal.FlatgKK;
    double* resultFlatgEtaPi= &theParamVal.FlatgEtaPi;
    
    double currentMass=par[resultCount];
    resultCount++;
    (*resultFlatMa980)=currentMass;
    
    double currentFlatgKK=par[resultCount];
    resultCount++;
    (*resultFlatgKK)=currentFlatgKK;
    
    double currentFlatgEtaPi=par[resultCount];
    resultCount++;
    (*resultFlatgEtaPi)=currentFlatgEtaPi;
  }
  else if (key==paramEnumKpKmPiGam::name(paramEnumKpKmPiGam::phaseSpace)){
    double* resultPhaseSpace= &theParamVal.phaseSpace;
    double currentPhaseSpace=par[resultCount];
    resultCount++;
    (*resultPhaseSpace)=currentPhaseSpace;
  }
  else { Alert << "Key: " << key << " not supported for setting up the Flatte mass parameters!!!" << endmsg;
    exit(1);
  }


  return resultCount;
}

void FitParamsKpKmPiGam::setMnUsrParamsDec(MnUserParameters& upar, paramKpKmPiGam& startVal,  paramKpKmPiGam& errVal, unsigned int index){
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
    std::string magStr=theJPCLS->name()+paramEnumKpKmPiGam::name(index)+"Mag";
    std::string phiStr=theJPCLS->name()+paramEnumKpKmPiGam::name(index)+"phi";

    double magVal=startPair.first;
    double phiVal=startPair.second;
    
    double magErr=errPair.first;
    double phiErr=errPair.second;
    
    double magMin=magVal-magErr;
    if (magMin<0.) magMin=0.;
    
    upar.Add(magStr, magVal, magErr, magMin, magVal+magErr);
    upar.Add(phiStr, phiVal, phiErr, -3.*M_PI, 3.*M_PI);

    counter++;
  }

}


void FitParamsKpKmPiGam::setMnUsrParamsMass(MnUserParameters& upar, paramKpKmPiGam& startVal,  paramKpKmPiGam& errVal, unsigned int index){
  pair<double, double> startParams;
  pair<double, double> errParams;
  startParams=massPair(startVal, index);
  errParams=massPair(errVal, index);

  //now fill the fitParameterMap
  std::string massStr=paramEnumKpKmPiGam::name(index)+"Mag";
  std::string widthStr=paramEnumKpKmPiGam::name(index)+"Width";
  
  double massVal=startParams.first;
  double widthVal=startParams.second;
  
  double massErr=errParams.first;
  double widthErr=errParams.second;

  double massMin=massVal-widthVal;
  if (massMin<0.) massMin=0.;

  double massMax=massVal+widthVal;

  double widthMin=0.; 
  double widthMax=2*widthVal;
 
  upar.Add(massStr, massVal, massErr, massMin, massMax);
  upar.Add(widthStr, widthVal, widthErr, widthMin, widthMax);
}


void FitParamsKpKmPiGam::setMnUsrParamsFlattea980Mass(MnUserParameters& upar, paramKpKmPiGam& startVal,  paramKpKmPiGam& errVal, std::string key){
 
  if (key==paramEnumKpKmPiGam::name(paramEnumKpKmPiGam::a980))
    {
      double start_a980M=startVal.FlatMa980;
      double start_FlatgKK=startVal.FlatgKK;
      double start_FlatgEtaPi=startVal.FlatgEtaPi;
      double err_a980M=errVal.FlatMa980;
      double err_FlatgKK=errVal.FlatgKK;
      double err_FlatgEtaPi=errVal.FlatgEtaPi;
      
      //now fill the fitParameterMap
      std::string massStr="a980Mass";
      std::string gKKStr="gKK";
      std::string gEtaPiStr="gEtaPi";
      
      double massMin=0.97;
      double massMax=1.2;
      upar.Add(massStr, start_a980M, err_a980M, massMin, massMax);
      
      upar.Add(gKKStr, start_FlatgKK, err_FlatgKK, 0., start_FlatgKK+3*err_FlatgKK);
      upar.Add(gEtaPiStr, start_FlatgEtaPi, err_FlatgEtaPi, 0., start_FlatgEtaPi+3*err_FlatgEtaPi);
    }
  
  else if ( key == paramEnumKpKmPiGam::name(paramEnumKpKmPiGam::phaseSpace)){ 
    double start_phaseSpace=startVal.phaseSpace;
    double err_phaseSpace=errVal.phaseSpace;
    std::string theStr="phaseSpace";
    
    double phaseSpaceMax=start_phaseSpace+3.*err_phaseSpace;
    if ( phaseSpaceMax > 1. ) phaseSpaceMax=1.;
    upar.Add(theStr, start_phaseSpace, err_phaseSpace, 0., phaseSpaceMax);
  }
  
  else { Alert << "Key: " << key << " not supported for setting up the Flatte mass parameters!!!" << endmsg;
    exit(1);
  } 
}
