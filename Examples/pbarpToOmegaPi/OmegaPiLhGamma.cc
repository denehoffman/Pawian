#include <getopt.h>
#include <fstream>
#include <string>


#include "Examples/pbarpToOmegaPi/OmegaPiLhGamma.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "PwaUtils/pbarpStates.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "ErrLogger/ErrLogger.hh"

#include <geneva/GConstrainedDoubleObject.hpp>
// #include <geneva/GConstrainedDoubleObjectCollection.hpp>
// #include <geneva/GDoubleGaussAdaptor.hpp>
// #include <geneva/GObjectExpectationChecksT.hpp>

OmegaPiLhGamma::OmegaPiLhGamma(boost::shared_ptr<const OmegaPiEventList> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0States> theStates) :
  AbsOmegaPiLh(theEvtList, theStates)
{
}

OmegaPiLhGamma::OmegaPiLhGamma(boost::shared_ptr<OmegaPiLhGamma> theOmegaPiLhGammaPtr):
  AbsOmegaPiLh(theOmegaPiLhGammaPtr)
{
}

OmegaPiLhGamma::~OmegaPiLhGamma()
{
}

double OmegaPiLhGamma::calcLogLh(const OmegaPiData::fitParamVal& theParamVal){

  double result=AbsOmegaPiLh::calcLogLh(theParamVal);
  if (_globalItCounter%100 == 0) printFitParams(std::cout, theParamVal);

//   if (_globalItCounter%1000 == 0){
//     std::ofstream theStream ("currentResult.dat");
//     dumpCurrentResult(theStream, theParamVal);
//     theStream.close();
//   }

  return result; 
}


double OmegaPiLhGamma::calcEvtIntensity(OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal){

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmSinglet=theParamVal.omegaProdSinglet;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet0=theParamVal.omegaProdTriplet0;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet1=theParamVal.omegaProdTriplet1;

  Spin lamgamma=-1;
  complex<double> singletAmpGM1=calcCoherentAmp(lamgamma,0, fitParmSinglet, theData);  
  complex<double> triplet0AmpGM1=calcCoherentAmp(lamgamma,0, fitParmTriplet0, theData);
  complex<double> tripletP1AmpGM1=calcCoherentAmp(lamgamma,1, fitParmTriplet1, theData);
  complex<double> tripletM1AmpGM1=calcCoherentAmp(lamgamma,-1, fitParmTriplet1, theData);

  lamgamma=1;
  complex<double> singletAmpGP1=calcCoherentAmp(lamgamma,0, fitParmSinglet, theData);  
  complex<double> triplet0AmpGP1=calcCoherentAmp(lamgamma,0, fitParmTriplet0, theData);
  complex<double> tripletP1AmpGP1=calcCoherentAmp(lamgamma,1, fitParmTriplet1, theData);
  complex<double> tripletM1AmpGP1=calcCoherentAmp(lamgamma,-1, fitParmTriplet1, theData);

  double result=norm(singletAmpGM1)+norm(triplet0AmpGM1)+norm(tripletP1AmpGM1)+norm(tripletM1AmpGM1);
  
  result+=norm(singletAmpGP1)+norm(triplet0AmpGP1)+norm(tripletP1AmpGP1)+norm(tripletM1AmpGP1);

  return result;  
}


complex<double> OmegaPiLhGamma::calcCoherentAmp(Spin lamgamma, Spin Minit, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm, OmegaPiData::OmPiEvtData* theData){

  complex<double> result(0.,0.);

  for (Spin lamomega=-1; lamomega<=1; lamomega++){

    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
    
    for ( it=fitParm.begin(); it!=fitParm.end(); ++it){
      boost::shared_ptr<const JPCLS> theJPCLS=it->first;
      pair<double, double> fitPair=it->second;     
      if (fabs(lamomega)>theJPCLS->J) continue;
      
      complex<double> omegaDecAmp=Clebsch(1,0,1,lamgamma,1, lamgamma)*conj(theData->Dfd[1][lamomega][lamgamma]);// Clebsch(1,lamgamma,0,0,1, lamgamma)=1
      complex<double> omegaProdAmp=calcOmegaProdPartAmp(Minit, lamomega, theJPCLS, fitPair, theData);
      
      result += omegaProdAmp*omegaDecAmp;
    }
  }
  return result;
}

complex<double>  OmegaPiLhGamma::calcOmegaProdPartAmp(Spin Minit, Spin lamomega, boost::shared_ptr<const JPCLS> theJPCLS, pair<double, double> fitVal, OmegaPiData::OmPiEvtData* theData){
 complex<double> result(0.,0.);

 if (fabs(lamomega)>theJPCLS->J) return result;
 
 double theMag=fitVal.first;
 double thePhi=fitVal.second;
 complex<double> expiphi(cos(thePhi), sin(thePhi));
 
 result=sqrt(2*theJPCLS->L+1)*sqrt(3)*theMag*expiphi*Clebsch(theJPCLS->L,0,1, lamomega,theJPCLS->J, lamomega)*conj(theData->Dfp[theJPCLS->J][Minit][lamomega]); //Clebsch(1,lamomega,0,0,1,lamomega)=1 
 
 if (Minit==-1 && theJPCLS->P==-1 && theJPCLS->C==-1){
   int theJ=theJPCLS->J;
   if( theJ!=0 && (theJ%2==0 ) ){ //J=even => negativ sign, due to CGs of the pbar p production
     result=-result;
   }
   
 }
 
 return result;
}


complex<double>  OmegaPiLhGamma::calcOmegaProdAmp(Spin Minit, Spin lamomega, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm, OmegaPiData::OmPiEvtData* theData){

  complex<double> result(0.,0.);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=fitParm.begin(); it!=fitParm.end(); ++it){

    boost::shared_ptr<const JPCLS> theJPCLS=it->first;
    pair<double, double> fitPair=it->second;
    result+=calcOmegaProdPartAmp(Minit, lamomega, theJPCLS, fitPair, theData);
  }
    
  return result;
}

complex<double> OmegaPiLhGamma::spinDensity(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal){

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmSinglet=theParamVal.omegaProdSinglet;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet0=theParamVal.omegaProdTriplet0;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet1=theParamVal.omegaProdTriplet1;

  complex<double> MsingletAmpGM1=calcOmegaProdAmp(0, M, fitParmSinglet, theData);  
  complex<double> Mtriplet0AmpGM1=calcOmegaProdAmp(0, M, fitParmTriplet0, theData);
  complex<double> MtripletP1AmpGM1=calcOmegaProdAmp(1, M, fitParmTriplet1, theData);
  complex<double> MtripletM1AmpGM1=calcOmegaProdAmp(-1, M, fitParmTriplet1, theData);

  complex<double> M_singletAmpGM1=calcOmegaProdAmp(0, M_, fitParmSinglet, theData);  
  complex<double> M_triplet0AmpGM1=calcOmegaProdAmp(0, M_, fitParmTriplet0, theData);
  complex<double> M_tripletP1AmpGM1=calcOmegaProdAmp(1, M_, fitParmTriplet1, theData);
  complex<double> M_tripletM1AmpGM1=calcOmegaProdAmp(-1, M_, fitParmTriplet1, theData);

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
      MsingletAmpGM1=calcOmegaProdAmp(0,M1, fitParmSinglet, theData);  
      Mtriplet0AmpGM1=calcOmegaProdAmp(0,M1, fitParmTriplet0, theData);
      MtripletP1AmpGM1=calcOmegaProdAmp(1,M1, fitParmTriplet1, theData);
      MtripletM1AmpGM1=calcOmegaProdAmp(-1,M1, fitParmTriplet1, theData);
      theNorm += (norm(MsingletAmpGM1)+norm(Mtriplet0AmpGM1)+norm(MtripletP1AmpGM1)+norm(MtripletM1AmpGM1));
    }
  }


  return (result/theNorm);

}

complex<double> OmegaPiLhGamma::spinDensityOmegaFrame(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal){

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


void OmegaPiLhGamma::getFitParamVal(OmegaPiData::fitParamVal& theParamVal, const std::vector<double>& par) const{
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaSinglet = _omegaPi0StatesPtr->jpclsSinglet();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0 = _omegaPi0StatesPtr->jpclsTriplet0();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1 = _omegaPi0StatesPtr->jpclsTriplet1();

  if (par.size()< JPCLSOmegaSinglet.size()*2+JPCLSOmegaTriplet0.size()*2+JPCLSOmegaTriplet1.size()*2-3) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\tJPCLSOmegaSinglet.size()+JPCLSOmegaTriplet0.size()+JPCLSOmegaTriplet1.size()-3=" << 
      JPCLSOmegaSinglet.size()*2+JPCLSOmegaTriplet0.size()*2+JPCLSOmegaTriplet1.size()*2-3 << endmsg;
    exit(1);
  }

  unsigned int counter=0;
  for ( itJPCLS=JPCLSOmegaSinglet.begin(); itJPCLS!=JPCLSOmegaSinglet.end(); ++itJPCLS){
    //now fill the fitParameterMap
    double mag=par[counter];
    counter++;
    double phi=0.;
    if (counter>1){ phi=par[counter];
    counter++;
    }
    std::pair <double,double> tmpParameter=make_pair(mag,phi);
    theParamVal.omegaProdSinglet[(*itJPCLS)]=tmpParameter; 
  }

  for ( itJPCLS=JPCLSOmegaTriplet0.begin(); itJPCLS!=JPCLSOmegaTriplet0.end(); ++itJPCLS){
    //now fill the fitParameterMap
    double mag=par[counter];
    counter++;
    double phi=0.;
    if (counter>JPCLSOmegaSinglet.size()*2){ phi=par[counter];
    counter++;
    }
    std::pair <double,double> tmpParameter=make_pair(mag,phi);
    theParamVal.omegaProdTriplet0[(*itJPCLS)]=tmpParameter; 
  }

  for ( itJPCLS=JPCLSOmegaTriplet1.begin(); itJPCLS!=JPCLSOmegaTriplet1.end(); ++itJPCLS){
    //now fill the fitParameterMap
    double mag=par[counter];
    counter++;
    double phi=0.;
    if (counter>JPCLSOmegaSinglet.size()*2+JPCLSOmegaTriplet0.size()*2-1){ phi=par[counter];
    counter++;
    }
    std::pair <double,double> tmpParameter=make_pair(mag,phi);
    theParamVal.omegaProdTriplet1[(*itJPCLS)]=tmpParameter; 
  }


}

void OmegaPiLhGamma::setGenevaFitParamVal( boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObjectCollection> theGbdc_ptr ){

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaSinglet= _omegaPi0StatesPtr->jpclsSinglet();

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  
  int counter=0;
  for ( itJPCLS=JPCLSOmegaSinglet.begin(); itJPCLS!=JPCLSOmegaSinglet.end(); ++itJPCLS){
    //now fill the fitParameterMap
    
    boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObject> gbd_ptr(new Gem::Geneva::GConstrainedDoubleObject(0., 1.) ); //JPCLS magnitude
    theGbdc_ptr->push_back(gbd_ptr);
    
    if (counter>0){ 
      boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObject>  gbd_ptr(new Gem::Geneva::GConstrainedDoubleObject(-M_PI, M_PI) ); //JPCLS phi
      theGbdc_ptr->push_back(gbd_ptr);
    }
    counter++; 
  }


  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0= _omegaPi0StatesPtr->jpclsTriplet0();
  counter=0;
  for ( itJPCLS=JPCLSOmegaTriplet0.begin(); itJPCLS!=JPCLSOmegaTriplet0.end(); ++itJPCLS){
    //now fill the fitParameterMap
    
    boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObject> gbd_ptr(new Gem::Geneva::GConstrainedDoubleObject(0., 1.) ); //JPCLS magnitude
    theGbdc_ptr->push_back(gbd_ptr);
    
    if (counter>0){ 
      boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObject>  gbd_ptr(new Gem::Geneva::GConstrainedDoubleObject(-M_PI, M_PI) ); //JPCLS phi
      theGbdc_ptr->push_back(gbd_ptr);
    }
    counter++; 
  }

 
 std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1= _omegaPi0StatesPtr->jpclsTriplet1();
  counter=0;
  for ( itJPCLS=JPCLSOmegaTriplet1.begin(); itJPCLS!=JPCLSOmegaTriplet1.end(); ++itJPCLS){
    //now fill the fitParameterMap
    
    boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObject> gbd_ptr(new Gem::Geneva::GConstrainedDoubleObject(0., 1.) ); //JPCLS magnitude
    theGbdc_ptr->push_back(gbd_ptr);
    
    if (counter>0){ 
      boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObject>  gbd_ptr(new Gem::Geneva::GConstrainedDoubleObject(-M_PI, M_PI) ); //JPCLS phi
      theGbdc_ptr->push_back(gbd_ptr);
    }
    counter++; 
  }


}


void OmegaPiLhGamma::setMnUsrParams(MnUserParameters& upar){

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaSinglet=_omegaPi0StatesPtr->jpclsSinglet();

  int counter=0;
  for ( itJPCLS=JPCLSOmegaSinglet.begin(); itJPCLS!=JPCLSOmegaSinglet.end(); ++itJPCLS){
    //now fill the fitParameterMap
    std::string magStr=(*itJPCLS)->name()+"S"+"mag";
    std::string phiStr=(*itJPCLS)->name()+"S"+"phi";

    upar.Add(magStr, 0.5, .1, 0., 2.);
    if (counter>0) upar.Add(phiStr, 0., .1, -3.*M_PI, 3.*M_PI);
    counter++; 
  }

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0=_omegaPi0StatesPtr->jpclsTriplet0();
  counter=0;
  for ( itJPCLS=JPCLSOmegaTriplet0.begin(); itJPCLS!=JPCLSOmegaTriplet0.end(); ++itJPCLS){
    //now fill the fitParameterMap
    std::string magStr=(*itJPCLS)->name()+"T0"+"mag";
    std::string phiStr=(*itJPCLS)->name()+"T0"+"phi";

    upar.Add(magStr, 0.5, .1, 0., 2.);
    if (counter>0) upar.Add(phiStr, 0., .1, -3.*M_PI, 3.*M_PI);
    counter++; 
  }

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1=_omegaPi0StatesPtr->jpclsTriplet1();
  counter=0;
  for ( itJPCLS=JPCLSOmegaTriplet1.begin(); itJPCLS!=JPCLSOmegaTriplet1.end(); ++itJPCLS){
    //now fill the fitParameterMap
    std::string magStr=(*itJPCLS)->name()+"T1"+"mag";
    std::string phiStr=(*itJPCLS)->name()+"T1"+"phi";

    upar.Add(magStr, 0.5, .1, 0., 2.);
    if (counter>0) upar.Add(phiStr, 0., .1, -3.*M_PI, 3.*M_PI);
    counter++; 
  }

}

void OmegaPiLhGamma::setMnUsrParams(MnUserParameters& upar, OmegaPiData::fitParamVal &finalFitParm){

   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmSinglet=finalFitParm.omegaProdSinglet;
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet0=finalFitParm.omegaProdTriplet0;
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet1=finalFitParm.omegaProdTriplet1;
   FillUpars(upar,fitParmSinglet,"S");
   FillUpars(upar,fitParmTriplet0,"T0");
   FillUpars(upar,fitParmTriplet1,"T1");
}

void OmegaPiLhGamma::FillUpars(MnUserParameters& upar, 
			       std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParmS,
			       const string &theSuffix)
{   
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
    bool bFirst=true;
    
    for ( it=fitParmS.begin(); it!=fitParmS.end(); ++it)
      {
	boost::shared_ptr<const JPCLS> theJPCLS=it->first;
	string strName = theJPCLS->name()+theSuffix;
	double theMag=it->second.first;
	double thePhi=it->second.second;
	
	//now fill the fitParameterMap
	std::string magStr=strName+"mag";
	std::string phiStr=strName+"phi";
	
	upar.Add(magStr, theMag, .1, 0., 2.);
	if (!bFirst) upar.Add(phiStr, thePhi, .1, -3.*M_PI, 3.*M_PI);
	if(bFirst) bFirst=false;
      } 
}

void  OmegaPiLhGamma::setMnUsrParams(MnUserParameters& upar, minuitStartParam &theStartParam)
{
  FillUpars(upar,theStartParam,_omegaPi0StatesPtr->jpclsSinglet(),"S");
  FillUpars(upar,theStartParam,_omegaPi0StatesPtr->jpclsTriplet0(),"T0");
  FillUpars(upar,theStartParam,_omegaPi0StatesPtr->jpclsTriplet1(),"T1");
}

//This Method fills the Minuit User Parameters with parameters for given state and if supllied initializes
//the parameter with start parameters given by user in theStartParam map.
void OmegaPiLhGamma::FillUpars(MnUserParameters& upar, 
		      minuitStartParam &theStartParam,
		      const std::vector< boost::shared_ptr<const JPCLS> > &theJPCLS,
		      const string &theSuffix
		      )
{   //now fill the fitParameterMap
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    bool bFirst=true;
    
    for ( itJPCLS=theJPCLS.begin(); itJPCLS!=theJPCLS.end(); ++itJPCLS)
    {
      string strName = (*itJPCLS)->name()+theSuffix;
      std::string magStr=strName+"mag";
      std::string phiStr=strName+"phi";
      ParameterMap::iterator it;

      it = theStartParam.getParamMap().find(strName);
      if (it != theStartParam.getParamMap().end()) 
      {
	upar.Add(magStr, it->second[0], .1, 0., 1.);
	if (!bFirst) upar.Add(phiStr, it->second[1], .1, -3.*M_PI, 3.*M_PI);
	DebugMsg << "\nUsing user start Parameter for the state " << strName << " with mag= " << it->second[0] 
	     << " and phi=" << it->second[1] << "\n"; 
      }
      else
      {
	upar.Add(magStr, 0.5, .1, 0., 1.);
	if (!bFirst) upar.Add(phiStr, 0., .1, -3.*M_PI, 3.*M_PI);
      }
      if (bFirst) bFirst=false;
    }
}

void OmegaPiLhGamma::dumpCurrentResult(std::ostream& os, const OmegaPiData::fitParamVal& fitParmVal) const{
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmSinglet=fitParmVal.omegaProdSinglet;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet0=fitParmVal.omegaProdTriplet0;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet1=fitParmVal.omegaProdTriplet1;

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::const_iterator iter;

  // singlet
  for ( iter=fitParmSinglet.begin(); iter!=fitParmSinglet.end(); ++iter){
    boost::shared_ptr<const JPCLS> theJPCLS=iter->first;
    double theMag=iter->second.first;
    double thePhi=iter->second.second;
    std::string strName = theJPCLS->name()+"S";
    os << strName << "\t" << theMag << "\t" << thePhi << "\n"; 
  }

  // triplet0
  for ( iter=fitParmTriplet0.begin(); iter!=fitParmTriplet0.end(); ++iter){
    boost::shared_ptr<const JPCLS> theJPCLS=iter->first;
    double theMag=iter->second.first;
    double thePhi=iter->second.second;
    std::string strName = theJPCLS->name()+"T0";
    os << strName << "\t" << theMag << "\t" << thePhi << "\n"; 
  }


  // triplet1
  for ( iter=fitParmTriplet1.begin(); iter!=fitParmTriplet1.end(); ++iter){
    boost::shared_ptr<const JPCLS> theJPCLS=iter->first;
    double theMag=iter->second.first;
    double thePhi=iter->second.second;
    std::string strName = theJPCLS->name()+"T1";
    os << strName << "\t" << theMag << "\t" << thePhi << "\n"; 
  }

}



void OmegaPiLhGamma::printFitParams(std::ostream& os, const OmegaPiData::fitParamVal& fitParmVal){

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmSinglet=fitParmVal.omegaProdSinglet;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet0=fitParmVal.omegaProdTriplet0;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet1=fitParmVal.omegaProdTriplet1;

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::const_iterator iter;

  os << "***fit parameter singlet states*** " << "\n";  
  for ( iter=fitParmSinglet.begin(); iter!=fitParmSinglet.end(); ++iter){
    os << iter->first->name()<< "\t";
    std::pair<double, double> tmpParam= iter->second;
    os <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << "\n";
  }
  
  os << "***fit parameter triplet m=0 states*** " << "\n";  
  for ( iter=fitParmTriplet0.begin(); iter!=fitParmTriplet0.end(); ++iter){
    os << iter->first->name()<< "\t";
    std::pair<double, double> tmpParam= iter->second;
    os <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << "\n";
  }
  
  os << "***fit parameter triplet m=1 states*** " << "\n";  
  for ( iter=fitParmTriplet1.begin(); iter!=fitParmTriplet1.end(); ++iter){
    os << iter->first->name()<< "\t";
    std::pair<double, double> tmpParam= iter->second;
    os <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << "\n";
  }
  
  
}

void OmegaPiLhGamma::print(std::ostream& os) const{
  os << "OmegaPiLhGamma\n";
}
