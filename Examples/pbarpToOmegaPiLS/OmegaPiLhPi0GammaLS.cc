#include <getopt.h>
#include <fstream>
#include <string>


#include "Examples/pbarpToOmegaPiLS/OmegaPiLhPi0GammaLS.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"
#include "PwaUtils/pbarpStatesLS.hh"
#include "Examples/pbarpToOmegaPiLS/pbarpToOmegaPi0StatesLS.hh"
#include "ErrLogger/ErrLogger.hh"

#include <geneva/GConstrainedDoubleObject.hpp>


OmegaPiLhPi0GammaLS::OmegaPiLhPi0GammaLS(boost::shared_ptr<const AbsOmegaPiEventListLS> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0StatesLS> theStates) :
  AbsOmegaPiLhLS(theEvtList, theStates)
{
}

OmegaPiLhPi0GammaLS::OmegaPiLhPi0GammaLS(boost::shared_ptr<OmegaPiLhPi0GammaLS> theOmegaPiLhPi0GammaLSPtr):
  AbsOmegaPiLhLS(theOmegaPiLhPi0GammaLSPtr)
{
}

OmegaPiLhPi0GammaLS::~OmegaPiLhPi0GammaLS()
{
}

double OmegaPiLhPi0GammaLS::calcLogLh(const OmegaPiDataLS::fitParamVal& theParamVal){

  double result=AbsOmegaPiLhLS::calcLogLh(theParamVal);
  if (_globalItCounter%10 == 0) printFitParams(std::cout, theParamVal);

//   if (_globalItCounter%1000 == 0){
//     std::ofstream theStream ("currentResult.dat");
//     dumpCurrentResult(theStream, theParamVal);
//     theStream.close();
//   }

  return result; 
}


double OmegaPiLhPi0GammaLS::calcEvtIntensity(OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal){

  Spin lamgamma=-1;
  complex<double> singletAmpGM1=calcCoherentAmp(lamgamma,0, theParamVal, _singlet_JPCLS_States, theData);  
  complex<double> triplet0AmpGM1=calcCoherentAmp(lamgamma,0, theParamVal, _triplet0_JPCLS_States, theData);
  complex<double> tripletP1AmpGM1=calcCoherentAmp(lamgamma,1, theParamVal, _tripletp1_JPCLS_States, theData);
  complex<double> tripletM1AmpGM1=calcCoherentAmp(lamgamma,-1, theParamVal, _tripletm1_JPCLS_States, theData);

  lamgamma=1;
  complex<double> singletAmpGP1=calcCoherentAmp(lamgamma,0, theParamVal, _singlet_JPCLS_States, theData);  
  complex<double> triplet0AmpGP1=calcCoherentAmp(lamgamma,0, theParamVal, _triplet0_JPCLS_States, theData);
  complex<double> tripletP1AmpGP1=calcCoherentAmp(lamgamma,1, theParamVal, _tripletp1_JPCLS_States, theData);
  complex<double> tripletM1AmpGP1=calcCoherentAmp(lamgamma,-1, theParamVal, _tripletm1_JPCLS_States, theData);

  double result=2.*norm(singletAmpGM1)+2.*norm(triplet0AmpGM1)+norm(tripletP1AmpGM1)+norm(tripletM1AmpGM1);
  
  result+=2.*norm(singletAmpGP1)+2.*norm(triplet0AmpGP1)+norm(tripletP1AmpGP1)+norm(tripletM1AmpGP1);

  return result;  
}


complex<double> OmegaPiLhPi0GammaLS::calcCoherentAmp(Spin lamgamma, Spin Minit, const OmegaPiDataLS::fitParamVal& theParamVal, std::vector< boost::shared_ptr<const JPCLSls> >& theJPCLSlsStates, OmegaPiDataLS::OmPiEvtDataLS* theData){

  complex<double> result(0.,0.);

  for (Spin lamomega=-1; lamomega<=1; lamomega++){

    complex<double> omegaDecAmp=Clebsch(1,0,1,lamgamma,1, lamgamma)*conj(theData->Dfd[1][lamomega][lamgamma]);// Clebsch(1,lamgamma,0,0,1, lamgamma)=1
    

    std::vector< boost::shared_ptr<const JPCLSls> >::const_iterator it;
     for ( it=theJPCLSlsStates.begin(); it!=theJPCLSlsStates.end(); ++it){

//        pair<double, double> fitPair= theParamVal.lsParam[(*it)];

       std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess >::const_iterator itmap;
       itmap = theParamVal.lsParam.find((*it));
       pair<double, double> fitPair= (*itmap).second; 
      
       if (fabs(lamomega)> (*it)->J ) continue;

       complex<double> omegaProdAmp=calcOmegaProdPartAmp(Minit, lamomega, (*it), fitPair, theData);
      
       result += omegaProdAmp*omegaDecAmp;
     }

  }

  return result;
  }

complex<double>  OmegaPiLhPi0GammaLS::calcOmegaProdPartAmp(Spin Minit, Spin lamomega, boost::shared_ptr<const JPCLSls> theJPCLSls, pair<double, double> fitVal, OmegaPiDataLS::OmPiEvtDataLS* theData){

 complex<double> result(0.,0.);

 if (fabs(lamomega)>theJPCLSls->J) return result;
 
 double theMag=fitVal.first;
 double thePhi=fitVal.second;
 complex<double> expiphi(cos(thePhi), sin(thePhi));
 
 result=sqrt(0.5)*theJPCLSls->preFactor*sqrt(2.*theJPCLSls->l+1)*Clebsch(theJPCLSls->l, 0, 1, lamomega, theJPCLSls->J, lamomega)*theMag*expiphi*conj(theData->Dfp[theJPCLSls->J][Minit][lamomega]);  //Clebsch(1, lamomega, 0, 0, theJPCLSls->s, lamomega)=1

 return result;
}


complex<double>  OmegaPiLhPi0GammaLS::calcOmegaProdAmp(Spin Minit, Spin lamomega, const OmegaPiDataLS::fitParamVal& theParamVal, std::vector< boost::shared_ptr<const JPCLSls> >& theJPCLSlsStates, OmegaPiDataLS::OmPiEvtDataLS* theData){

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

complex<double> OmegaPiLhPi0GammaLS::spinDensity(Spin M, Spin M_, OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal){

 
  complex<double> MsingletAmpGM1=calcOmegaProdAmp(0, M, theParamVal, _singlet_JPCLS_States, theData);  
  complex<double> Mtriplet0AmpGM1=calcOmegaProdAmp(0, M, theParamVal, _triplet0_JPCLS_States, theData);
  complex<double> MtripletP1AmpGM1=calcOmegaProdAmp(1, M, theParamVal, _tripletp1_JPCLS_States, theData);
  complex<double> MtripletM1AmpGM1=calcOmegaProdAmp(-1, M, theParamVal, _tripletm1_JPCLS_States, theData);

  complex<double> M_singletAmpGM1=calcOmegaProdAmp(0, M_, theParamVal, _singlet_JPCLS_States, theData);  
  complex<double> M_triplet0AmpGM1=calcOmegaProdAmp(0, M_, theParamVal, _triplet0_JPCLS_States, theData);
  complex<double> M_tripletP1AmpGM1=calcOmegaProdAmp(1, M_, theParamVal, _tripletp1_JPCLS_States, theData);
  complex<double> M_tripletM1AmpGM1=calcOmegaProdAmp(-1, M_, theParamVal, _tripletm1_JPCLS_States, theData);

  complex<double> result(0.,0.);

  result=2.*MsingletAmpGM1*conj(M_singletAmpGM1)
    +2.*Mtriplet0AmpGM1*conj(M_triplet0AmpGM1)
    +MtripletP1AmpGM1*conj(M_tripletP1AmpGM1)
    +MtripletM1AmpGM1*conj(M_tripletM1AmpGM1);
  
    double theNorm = 2.*norm(MsingletAmpGM1)+2.*norm(Mtriplet0AmpGM1)+norm(MtripletP1AmpGM1)+norm(MtripletM1AmpGM1);


    for (Spin M1=-1; M1<=1; M1++)
  {
    if(M1!=M)
    {
      MsingletAmpGM1=calcOmegaProdAmp(0,M1, theParamVal, _singlet_JPCLS_States, theData);  
      Mtriplet0AmpGM1=calcOmegaProdAmp(0,M1, theParamVal, _triplet0_JPCLS_States, theData);
      MtripletP1AmpGM1=calcOmegaProdAmp(1,M1, theParamVal, _tripletp1_JPCLS_States, theData);
      MtripletM1AmpGM1=calcOmegaProdAmp(-1,M1, theParamVal, _tripletm1_JPCLS_States, theData);
      theNorm += (2.*norm(MsingletAmpGM1)+2.*norm(Mtriplet0AmpGM1)+norm(MtripletP1AmpGM1)+norm(MtripletM1AmpGM1));
    }
  }


 return (result/theNorm);

}

complex<double> OmegaPiLhPi0GammaLS::spinDensityOmegaFrame(Spin M, Spin M_, OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal){

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


void OmegaPiLhPi0GammaLS::getFitParamVal(OmegaPiDataLS::fitParamVal& theParamVal, const std::vector<double>& par) const{
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

void OmegaPiLhPi0GammaLS::setGenevaFitParamVal( boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObjectCollection> theGbdc_ptr ){

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


void OmegaPiLhPi0GammaLS::setMnUsrParams(MnUserParameters& upar){

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

void OmegaPiLhPi0GammaLS::setMnUsrParams(MnUserParameters& upar, OmegaPiDataLS::fitParamVal &fitParam){


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



void  OmegaPiLhPi0GammaLS::setMnUsrParams(MnUserParameters& upar, MinuitStartParamLS &theStartParam)
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


void OmegaPiLhPi0GammaLS::dumpCurrentResult(std::ostream& os, const OmegaPiDataLS::fitParamVal& fitParmVal) const{
  std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess > lsFitParam=fitParmVal.lsParam;

  std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess >::const_iterator iter;
  
  for ( iter=lsFitParam.begin(); iter != lsFitParam.end(); ++iter){
    os << iter->first->name()<< "\t";
    std::pair<double, double> tmpParam= iter->second;
    os <<"\t" << tmpParam.first <<"\t" << tmpParam.second  << "\n";
  }

}



void OmegaPiLhPi0GammaLS::printFitParams(std::ostream& os, const OmegaPiDataLS::fitParamVal& fitParmVal){

  std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess > lsFitParam=fitParmVal.lsParam;

  std::map< boost::shared_ptr<const JPCLSls>, pair<double, double>, pawian::Collection::SharedPtrLess >::const_iterator iter;

  os << "\n***fit parameter *** " << "\n";  
  for ( iter=lsFitParam.begin(); iter != lsFitParam.end(); ++iter){
    os << iter->first->name()<< "\t";
    std::pair<double, double> tmpParam= iter->second;
    os <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << "\n";
  }
}

void OmegaPiLhPi0GammaLS::print(std::ostream& os) const{
  os << "OmegaPiLhPi0GammaLS\n";
}
