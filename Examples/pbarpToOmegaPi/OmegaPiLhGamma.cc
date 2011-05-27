#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarpToOmegaPi/OmegaPiLhGamma.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "PwaUtils/pbarpStates.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "ErrLogger/ErrLogger.hh"

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



void OmegaPiLhGamma::print(std::ostream& os) const{
  os << "OmegaPiLhGamma::print\n";
}
