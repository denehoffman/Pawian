#include "Examples/pbarpToOmegaPi/spindensity.hh"

#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "PwaUtils/pbarpStates.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "ErrLogger/ErrLogger.hh"

spinDensity::spinDensity()
{
}


spinDensity::~spinDensity()
{
}


complex<double> spinDensity::calcSpinDensityIncoherent(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal)
{
//   double thephi=theData->omegaHeliCm4Vec.Phi();spindensity.cc
//   double thetheta=theData->omegaHeliCm4Vec.Theta();
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmSinglet=theParamVal.omegaProdSinglet;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet0=theParamVal.omegaProdTriplet0;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet1=theParamVal.omegaProdTriplet1;

  complex<double> MsingletAmpGM1=calcIncoherentAmp(0, M, fitParmSinglet, theData);  
  complex<double> Mtriplet0AmpGM1=calcIncoherentAmp(0, M, fitParmTriplet0, theData);
  complex<double> MtripletP1AmpGM1=calcIncoherentAmp(1, M, fitParmTriplet1, theData);
  complex<double> MtripletM1AmpGM1=calcIncoherentAmp(-1, M, fitParmTriplet1, theData);

  complex<double> M_singletAmpGM1=calcIncoherentAmp(0, M_, fitParmSinglet, theData);  
  complex<double> M_triplet0AmpGM1=calcIncoherentAmp(0, M_, fitParmTriplet0, theData);
  complex<double> M_tripletP1AmpGM1=calcIncoherentAmp(1, M_, fitParmTriplet1, theData);
  complex<double> M_tripletM1AmpGM1=calcIncoherentAmp(-1, M_, fitParmTriplet1, theData);

   complex<double> result;
   result=MsingletAmpGM1*conj(M_singletAmpGM1)
     +Mtriplet0AmpGM1*conj(M_triplet0AmpGM1)
     +MtripletP1AmpGM1*conj(M_tripletP1AmpGM1)
     +MtripletM1AmpGM1*conj(M_tripletM1AmpGM1);
   


    double N = norm(MsingletAmpGM1)+norm(Mtriplet0AmpGM1)+norm(MtripletP1AmpGM1)+norm(MtripletM1AmpGM1);


    for (Spin M1=-1; M1<=1; M1++)
  {
    if(M1!=M)
    {
      MsingletAmpGM1=calcIncoherentAmp(0,M1, fitParmSinglet, theData);  
      Mtriplet0AmpGM1=calcIncoherentAmp(0,M1, fitParmTriplet0, theData);
      MtripletP1AmpGM1=calcIncoherentAmp(1,M1, fitParmTriplet1, theData);
      MtripletM1AmpGM1=calcIncoherentAmp(-1,M1, fitParmTriplet1, theData);
      N += (norm(MsingletAmpGM1)+norm(Mtriplet0AmpGM1)+norm(MtripletP1AmpGM1)+norm(MtripletM1AmpGM1));
    }
  }


  return (result/N);  
}

complex<double> spinDensity::calcSpinDensityOmegaFrame(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal){

  complex<double> rhoAdair=spinDensity::calcSpinDensityIncoherent(M, M_, theData, theParamVal);
  double thetaOmegaCms=theData->omegaHeliCm4Vec.Theta();
  complex<double> result(0.,0.);

  for (Spin i=-1; i<=1; i++){
    for (Spin j=-1; j<=1; j++){
      complex<double> rhoAdair=calcSpinDensityIncoherent(i, j, theData, theParamVal);
      result+=Wigner_d(1,M,i,-thetaOmegaCms)*rhoAdair*Wigner_d(1,j,M_,thetaOmegaCms);
    }
  }

  return result;
}



complex<double> spinDensity::calcIncoherentAmp(Spin Minit, Spin lamomega, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm, OmegaPiData::OmPiEvtData* theData)
{

  complex<double> result(0.,0.);
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=fitParm.begin(); it!=fitParm.end(); ++it)
  {
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;
     if (fabs(lamomega)>theJPCLS->J) continue;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    complex<double> tmpResult=sqrt(2*theJPCLS->L+1)*theMag*expiphi*Clebsch(theJPCLS->L,0,1, lamomega,theJPCLS->J, lamomega)*conj(theData->Dfp[theJPCLS->J][Minit][lamomega]); //Clebsch(1,lamomega,0,0,1,lamomega)=1

     
    if (Minit==-1 && theJPCLS->P==-1 && theJPCLS->C==-1){
      int theJ=theJPCLS->J;
      if( theJ!=0 && (theJ%2==0 ) ){ //J=even => negativ sign, due to CGs of the pbar p production
        tmpResult=-tmpResult;
      }
    }

    result+=tmpResult;
    
  } 

  return result;
}

