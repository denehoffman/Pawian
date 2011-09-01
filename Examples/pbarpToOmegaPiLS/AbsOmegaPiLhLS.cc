#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiLhLS.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"
#include "PwaUtils/pbarpStatesLS.hh"
#include "Examples/pbarpToOmegaPiLS/pbarpToOmegaPi0StatesLS.hh"
#include "ErrLogger/ErrLogger.hh"
#include <geneva/GConstrainedDoubleObject.hpp>



AbsOmegaPiLhLS::AbsOmegaPiLhLS(boost::shared_ptr<const AbsOmegaPiEventListLS> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0StatesLS> theStates) :
  _globalItCounter(0),
  _omegaPiEventListPtr(theEvtList),
  _omegaPi0StatesPtr(theStates)
{
  _evtDataVec=_omegaPiEventListPtr->getDataVecs();
  _evtMCVec=_omegaPiEventListPtr->getMcVecs();

  _allJPCLSlsStates=theStates->all_JPCLSls_States();
  _singlet_JPCLS_States=theStates->singlet_JPCLSls_States();
  _triplet0_JPCLS_States=theStates->triplet0_JPCLSls_States();
  _tripletp1_JPCLS_States=theStates->tripletp1_JPCLSls_States();
  _tripletm1_JPCLS_States=theStates->tripletm1_JPCLSls_States();
}

AbsOmegaPiLhLS::AbsOmegaPiLhLS(boost::shared_ptr<AbsOmegaPiLhLS> theAbsOmegaPiLhLSPtr):
  _omegaPiEventListPtr(theAbsOmegaPiLhLSPtr->getEventList()),
  _omegaPi0StatesPtr(theAbsOmegaPiLhLSPtr->omegaPi0States())
{
  _evtDataVec=_omegaPiEventListPtr->getDataVecs();
  _evtMCVec=_omegaPiEventListPtr->getMcVecs();

  _allJPCLSlsStates=_omegaPi0StatesPtr->all_JPCLSls_States();
  _singlet_JPCLS_States=_omegaPi0StatesPtr->singlet_JPCLSls_States();
  _triplet0_JPCLS_States= _omegaPi0StatesPtr->triplet0_JPCLSls_States();
  _tripletp1_JPCLS_States= _omegaPi0StatesPtr->tripletp1_JPCLSls_States();
  _tripletm1_JPCLS_States= _omegaPi0StatesPtr->tripletm1_JPCLSls_States();
}

AbsOmegaPiLhLS::~AbsOmegaPiLhLS()
{
}

double AbsOmegaPiLhLS::calcLogLh(const OmegaPiDataLS::fitParamVal& theParamVal){

  _globalItCounter++;
 
  double logLH=0.;
  double logLH_data=0.;
  double weightSum=0.;

  std::vector<OmegaPiDataLS::OmPiEvtDataLS*>::iterator iterd;
  for (iterd=_evtDataVec.begin(); iterd!=_evtDataVec.end(); ++iterd){
    double intensity=calcEvtIntensity((*iterd), theParamVal);
    if (intensity>0.) logLH_data+= ((*iterd)->eventWeight) * log10(intensity);
    weightSum+= (*iterd)->eventWeight;
  } 

  double LH_mc=0.;
  
  std::vector<OmegaPiDataLS::OmPiEvtDataLS*>::iterator iterm;
  for (iterm=_evtMCVec.begin(); iterm!=_evtMCVec.end(); ++iterm){
           double intensity=calcEvtIntensity((*iterm), theParamVal);
           LH_mc+=intensity;
         }

  double logLH_mc_Norm=0.;
  if (LH_mc>0.) logLH_mc_Norm=log10(LH_mc/_evtMCVec.size());

  logLH= weightSum *(LH_mc/_evtMCVec.size()-1.)*(LH_mc/_evtMCVec.size()-1.)
    -2.*logLH_data
    +2.*weightSum*logLH_mc_Norm;

  Info << "current LH = " << logLH << endmsg;

 return logLH;

}

void AbsOmegaPiLhLS::getFitParamVal(OmegaPiDataLS::fitParamVal& theParamVal, const std::vector<double>& par) const{
  std::vector< boost::shared_ptr<const JPCLSls> >::const_iterator itJPCLSls;


  if (par.size()< _allJPCLSlsStates.size()*2-2) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t_allJPCLSlsStates.size()*2-2=" << 
      _allJPCLSlsStates.size()*2-2 << endmsg;
    exit(1);
  }

  unsigned int counter=0;
  for ( itJPCLSls=_allJPCLSlsStates.begin(); itJPCLSls!=_allJPCLSlsStates.end(); ++itJPCLSls){

    bool fillPhi=true;
    if( *(*itJPCLSls)== *(*_singlet_JPCLS_States.begin()) || *(*itJPCLSls)== *(*_triplet0_JPCLS_States.begin()) ) fillPhi=false;

    //now fill the fitParameterMap
    double mag=par[counter];
    counter++;
    double phi=0.;
    if (fillPhi){ phi=par[counter];
    counter++;
    }
    std::pair <double,double> tmpParameter=make_pair(mag,phi);
    theParamVal.lsParam[(*itJPCLSls)]=tmpParameter; 
  }

}

void AbsOmegaPiLhLS::setGenevaFitParamVal( boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObjectCollection> theGbdc_ptr ){

  std::vector< boost::shared_ptr<const JPCLSls> >::const_iterator itJPCLSls;

  for ( itJPCLSls=_allJPCLSlsStates.begin(); itJPCLSls!=_allJPCLSlsStates.end(); ++itJPCLSls){

    bool fixPhi=false;

    if( *(*itJPCLSls)== *(*_singlet_JPCLS_States.begin()) || *(*itJPCLSls)== *(*_triplet0_JPCLS_States.begin()) ) fixPhi=true;

      boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObject> gbd_ptr(new Gem::Geneva::GConstrainedDoubleObject(0., 1.) ); //JPCLSls magnitude
     theGbdc_ptr->push_back(gbd_ptr);
     if (!fixPhi){
      boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObject>  gbd_ptr(new Gem::Geneva::GConstrainedDoubleObject(-M_PI, M_PI) ); //JPCLS phi
      theGbdc_ptr->push_back(gbd_ptr);
     }

  }

}

void AbsOmegaPiLhLS::setMnUsrParams(MnUserParameters& upar){

  std::vector< boost::shared_ptr<const JPCLSls> >::const_iterator itJPCLSls;

  for ( itJPCLSls=_allJPCLSlsStates.begin(); itJPCLSls!=_allJPCLSlsStates.end(); ++itJPCLSls){

    bool fixPhi=false;

    if( *(*itJPCLSls)== *(*_singlet_JPCLS_States.begin()) || *(*itJPCLSls)== *(*_triplet0_JPCLS_States.begin()) ) fixPhi=true;

      string strName = (*itJPCLSls)->name();
      std::string magStr=strName+"mag";
      std::string phiStr=strName+"phi";

      upar.Add(magStr, 0.5, .1, 0., 1.);
      if (!fixPhi) upar.Add(phiStr, 0., .1, -3.*M_PI, 3.*M_PI);

  }

}


void AbsOmegaPiLhLS::setMnUsrParams(MnUserParameters& upar, OmegaPiDataLS::fitParamVal &fitParam){


  std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;

  for ( it=fitParam.lsParam.begin(); it!=fitParam.lsParam.end(); ++it){

    bool fixPhi=false;
    if( *(it->first) == *(*_singlet_JPCLS_States.begin()) || *(it->first)== *(*_triplet0_JPCLS_States.begin()) ) fixPhi=true;

    string strName = it->first->name();
    std::string magStr=strName+"mag";
    std::string phiStr=strName+"phi";

    double theMag=it->second.first;
    double thePhi=it->second.second;

    upar.Add(magStr, theMag, .1, 0., 1.);
    if (!fixPhi) upar.Add(phiStr, thePhi, .1, -3.*M_PI, 3.*M_PI);
  }

}


void  AbsOmegaPiLhLS::setMnUsrParams(MnUserParameters& upar, MinuitStartParamLS &theStartParam)
{

  std::vector< boost::shared_ptr<const JPCLSls> >::const_iterator itJPCLSls;

  for ( itJPCLSls=_allJPCLSlsStates.begin(); itJPCLSls!=_allJPCLSlsStates.end(); ++itJPCLSls){

    bool fixPhi=false;

    if( *(*itJPCLSls)== *(*_singlet_JPCLS_States.begin()) || *(*itJPCLSls)== *(*_triplet0_JPCLS_States.begin()) ) fixPhi=true;

      string strName = (*itJPCLSls)->name();
      std::string magStr=strName+"mag";
      std::string phiStr=strName+"phi";
      ParameterMap::iterator itMap;

      itMap = theStartParam.getParamMap().find(strName);

      if (itMap != theStartParam.getParamMap().end()) 
	{
	  upar.Add(magStr, itMap->second[0], .1, 0., 1.);
	  if (!fixPhi) upar.Add(phiStr, itMap->second[1], .1, -3.*M_PI, 3.*M_PI);
	  DebugMsg << "\nUsing user start Parameter for the state " << strName << " with mag= " << itMap->second[0] 
		   << " and phi=" << itMap->second[1] << "\n"; 
	}
      else
	{
	  upar.Add(magStr, 0.5, .1, 0., 1.);
	  if (!fixPhi) upar.Add(phiStr, 0., .1, -3.*M_PI, 3.*M_PI);
	}
  }
  
}


void AbsOmegaPiLhLS::fixMnUsrParams(MnUserParameters& upar, std::vector<std::string>& fixedParams){


  //first test wheather all fixParams exist
  std::vector<double> parmVecs=upar.Params();
  std::vector<std::string> uparNamesVec;

  for (unsigned int i=0; i<parmVecs.size(); ++i){
    uparNamesVec.push_back(upar.GetName(i));    
  }

  std::vector<std::string>::const_iterator it;
  for (it=fixedParams.begin(); it!=fixedParams.end(); ++it){
    std::string tmpName = (*it);
    std::vector<std::string>::const_iterator found=std::find( uparNamesVec.begin(), uparNamesVec.end(), tmpName );

    if ( (found == uparNamesVec.end()) ){
      Alert << "list for fixing MINUIT parameters contains wrong names!!!" << "\n"
	    << tmpName << " is not a MINUIT fit parameter name" << endmsg;
      exit(0);
    }
  }


  for (it=fixedParams.begin(); it != fixedParams.end(); ++it){
    DebugMsg << "Index(" << (*it) << "):\t" << upar.Index( (*it) ) << endmsg;
    upar.Fix( (*it) );
  }

}


void AbsOmegaPiLhLS::dumpCurrentResult(std::ostream& os, const OmegaPiDataLS::fitParamVal& fitParmVal) const{
  std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess > lsFitParam=fitParmVal.lsParam;

  std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess >::const_iterator iter;
  
  for ( iter=lsFitParam.begin(); iter != lsFitParam.end(); ++iter){
    os << iter->first->name()<< "\t";
    std::pair<double, double> tmpParam= iter->second;
    os <<"\t" << tmpParam.first <<"\t" << tmpParam.second  << "\n";
  }

}


complex<double>  AbsOmegaPiLhLS::calcOmegaProdPartAmp(Spin Minit, Spin lamomega, boost::shared_ptr<const JPCLSls> theJPCLSls, pair<double, double> fitVal, OmegaPiDataLS::OmPiEvtDataLS* theData){

 complex<double> result(0.,0.);

 if (fabs(lamomega)>theJPCLSls->J) return result;
 
 double theMag=fitVal.first;
 double thePhi=fitVal.second;
 complex<double> expiphi(cos(thePhi), sin(thePhi));
 
 result=theJPCLSls->preFactor*sqrt(2.*theJPCLSls->l+1)*Clebsch(theJPCLSls->l, 0, 1, lamomega, theJPCLSls->J, lamomega)*theMag*expiphi*theData->Dfp[theJPCLSls->J][Minit][lamomega];  //Clebsch(1, lamomega, 0, 0, theJPCLSls->s, lamomega)=1
 return result;
}


complex<double>  AbsOmegaPiLhLS::calcOmegaProdAmp(Spin Minit, Spin lamomega, const OmegaPiDataLS::fitParamVal& theParamVal, std::vector< boost::shared_ptr<const JPCLSls> >& theJPCLSlsStates, OmegaPiDataLS::OmPiEvtDataLS* theData){

  complex<double> result(0.,0.);

  std::vector< boost::shared_ptr<const JPCLSls> >::const_iterator it;
  for ( it=theJPCLSlsStates.begin(); it!=theJPCLSlsStates.end(); ++it){

    std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess >::const_iterator itmap;
    itmap = theParamVal.lsParam.find((*it));
    pair<double, double> fitPair= (*itmap).second; 
      
    if (fabs(lamomega)> (*it)->J ) continue;

    result+=calcOmegaProdPartAmp(Minit, lamomega, (*it), fitPair, theData);
  }
    
  return result;
}



complex<double> AbsOmegaPiLhLS::spinDensity(Spin M, Spin M_, OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal){

 
  complex<double> MsingletAmpGM1=calcOmegaProdAmp(0, M, theParamVal, _singlet_JPCLS_States, theData);  
  complex<double> Mtriplet0AmpGM1=calcOmegaProdAmp(0, M, theParamVal, _triplet0_JPCLS_States, theData);
  complex<double> MtripletP1AmpGM1=calcOmegaProdAmp(1, M, theParamVal, _tripletp1_JPCLS_States, theData);
  complex<double> MtripletM1AmpGM1=calcOmegaProdAmp(-1, M, theParamVal, _tripletm1_JPCLS_States, theData);

  complex<double> M_singletAmpGM1=calcOmegaProdAmp(0, M_, theParamVal, _singlet_JPCLS_States, theData);  
  complex<double> M_triplet0AmpGM1=calcOmegaProdAmp(0, M_, theParamVal, _triplet0_JPCLS_States, theData);
  complex<double> M_tripletP1AmpGM1=calcOmegaProdAmp(1, M_, theParamVal, _tripletp1_JPCLS_States, theData);
  complex<double> M_tripletM1AmpGM1=calcOmegaProdAmp(-1, M_, theParamVal, _tripletm1_JPCLS_States, theData);

  complex<double> result(0.,0.);

  result=MsingletAmpGM1*conj(M_singletAmpGM1)
    +Mtriplet0AmpGM1*conj(M_triplet0AmpGM1)
    +MtripletP1AmpGM1*conj(M_tripletP1AmpGM1)
    +MtripletM1AmpGM1*conj(M_tripletM1AmpGM1);
  
    double theNorm = norm(MsingletAmpGM1)+norm(Mtriplet0AmpGM1)+norm(MtripletP1AmpGM1)+norm(MtripletM1AmpGM1);


    for (Spin M1=-1; M1<=1; M1++)
  {
    if(M1!=M)
    {
      MsingletAmpGM1=calcOmegaProdAmp(0,M1, theParamVal, _singlet_JPCLS_States, theData);  
      Mtriplet0AmpGM1=calcOmegaProdAmp(0,M1, theParamVal, _triplet0_JPCLS_States, theData);
      MtripletP1AmpGM1=calcOmegaProdAmp(1,M1, theParamVal, _tripletp1_JPCLS_States, theData);
      MtripletM1AmpGM1=calcOmegaProdAmp(-1,M1, theParamVal, _tripletm1_JPCLS_States, theData);
      theNorm += (norm(MsingletAmpGM1)+norm(Mtriplet0AmpGM1)+norm(MtripletP1AmpGM1)+norm(MtripletM1AmpGM1));
    }
  }


 return (result/theNorm);

}

complex<double> AbsOmegaPiLhLS::spinDensityOmegaFrame(Spin M, Spin M_, OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal){

  double thetaOmegaCms=theData->omegaHeliCm4Vec.Theta();
  complex<double> result(0.,0.);

  for (Spin i=-1; i<=1; i++){
    for (Spin j=-1; j<=1; j++){
      complex<double> rhoAdair=spinDensity(i, j, theData, theParamVal);
      result+=Wigner_d(1,M,i,-thetaOmegaCms)*rhoAdair*Wigner_d(1,j,M_,thetaOmegaCms);
    }
  }

  return result;

}



void AbsOmegaPiLhLS::printFitParams(std::ostream& os, const OmegaPiDataLS::fitParamVal& fitParmVal){

  std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess > lsFitParam=fitParmVal.lsParam;

  std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess >::const_iterator iter;

  os << "\n***fit parameter *** " << "\n";  
  for ( iter=lsFitParam.begin(); iter != lsFitParam.end(); ++iter){
    os << iter->first->name()<< "\t";
    std::pair<double, double> tmpParam= iter->second;
    os <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << "\n";
  }
}


void AbsOmegaPiLhLS::print(std::ostream& os) const{
  os << "AbsOmegaPiLhLS::print\n";
}
