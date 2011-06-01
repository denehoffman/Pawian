#include <getopt.h>
#include <fstream>
#include <string>


#include "Examples/pbarpToOmegaPi/OmegaPiLhGammaBw.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "PwaUtils/pbarpStates.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "ErrLogger/ErrLogger.hh"

#include <geneva/GConstrainedDoubleObject.hpp>


OmegaPiLhGammaBw::OmegaPiLhGammaBw(boost::shared_ptr<const OmegaPiEventList> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0States> theStates) :
  OmegaPiLhGamma(theEvtList, theStates)
{
}

OmegaPiLhGammaBw::OmegaPiLhGammaBw(boost::shared_ptr<OmegaPiLhGammaBw> theOmegaPiLhGammaBwPtr):
  OmegaPiLhGamma(theOmegaPiLhGammaBwPtr)
{
}

OmegaPiLhGammaBw::~OmegaPiLhGammaBw()
{
}

double OmegaPiLhGammaBw::calcLogLh(const OmegaPiData::fitParamVal& theParamVal){

  double result=AbsOmegaPiLh::calcLogLh(theParamVal);
  if (_globalItCounter%100 == 0) printFitParams(std::cout, theParamVal);
//   if (_globalItCounter%1000 == 0){
//     std::ofstream theStream ("currentResult.dat");
//     dumpCurrentResult(theStream, theParamVal);
//     theStream.close();
//   }
  return result; 
}


double OmegaPiLhGammaBw::calcEvtIntensity(OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal){

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmSinglet=theParamVal.omegaProdSinglet;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet0=theParamVal.omegaProdTriplet0;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet1=theParamVal.omegaProdTriplet1;

  double omegaMassFit=theParamVal.omegaMass.first;
  double omegaWidthFit=theParamVal.omegaMass.second;


  Spin lamgamma=-1;

  complex<double> singletAmpGM1=calcCoherentAmp(lamgamma,0, fitParmSinglet, omegaMassFit, omegaWidthFit, theData); 
  complex<double> triplet0AmpGM1=calcCoherentAmp(lamgamma,0, fitParmTriplet0,omegaMassFit, omegaWidthFit, theData);
  complex<double> tripletP1AmpGM1=calcCoherentAmp(lamgamma,1, fitParmTriplet1,omegaMassFit, omegaWidthFit, theData);
  complex<double> tripletM1AmpGM1=calcCoherentAmp(lamgamma,-1, fitParmTriplet1,omegaMassFit, omegaWidthFit, theData);

  lamgamma=1;
  complex<double> singletAmpGP1=calcCoherentAmp(lamgamma,0, fitParmSinglet,omegaMassFit, omegaWidthFit, theData);
  complex<double> triplet0AmpGP1=calcCoherentAmp(lamgamma,0, fitParmTriplet0,omegaMassFit, omegaWidthFit, theData);
  complex<double> tripletP1AmpGP1=calcCoherentAmp(lamgamma,1, fitParmTriplet1,omegaMassFit, omegaWidthFit, theData);
  complex<double> tripletM1AmpGP1=calcCoherentAmp(lamgamma,-1, fitParmTriplet1,omegaMassFit, omegaWidthFit, theData);

  double result=norm(singletAmpGM1)+norm(triplet0AmpGM1)+norm(tripletP1AmpGM1)+norm(tripletM1AmpGM1);
  
  result+=norm(singletAmpGP1)+norm(triplet0AmpGP1)+norm(tripletP1AmpGP1)+norm(tripletM1AmpGP1);

  return result;  
}


complex<double> OmegaPiLhGammaBw::calcCoherentAmp(Spin lamgamma, Spin Minit, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm, double fitParmOmMass, double fitParmOmWidth, OmegaPiData::OmPiEvtData* theData){

  complex<double> result(0.,0.);

  for (Spin lamomega=-1; lamomega<=1; lamomega++){

    complex<double> omegaDecAmp1=Clebsch(1,0,1,lamgamma,1, lamgamma)*conj(theData->Dfd[1][lamomega][lamgamma]);// Clebsch(1,lamgamma,0,0,1, lamgamma)=1
    complex<double> omegaDecAmp2=Clebsch(1,0,1,lamgamma,1, lamgamma)*conj(theData->Dfd2[1][lamomega][lamgamma]);// Clebsch(1,lamgamma,0,0,1, lamgamma)=1    
    
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
    
    for ( it=fitParm.begin(); it!=fitParm.end(); ++it){
      boost::shared_ptr<const JPCLS> theJPCLS=it->first;
      pair<double, double> fitPair=it->second;     
      if (fabs(lamomega)>theJPCLS->J) continue;
      
      complex<double> omegaProdAmp1=calcOmegaProdPartAmp(Minit, lamomega, theJPCLS, fitPair, fitParmOmMass, fitParmOmWidth, theData, 1);
      complex<double> omegaProdAmp2=calcOmegaProdPartAmp(Minit, lamomega, theJPCLS, fitPair, fitParmOmMass, fitParmOmWidth, theData, 2);
      
      result += omegaProdAmp1*omegaDecAmp1+omegaProdAmp2*omegaDecAmp2;
    }
  }
  return result;
}

complex<double> OmegaPiLhGammaBw::calcOmegaProdAmp(Spin Minit, Spin lamomega, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm,  double fitParmOmMass, double fitParmOmWidth, OmegaPiData::OmPiEvtData* theData){

  complex<double> result(0.,0.);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=fitParm.begin(); it!=fitParm.end(); ++it){
    
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;
    pair<double, double> fitPair=it->second;
    result+=calcOmegaProdPartAmp(Minit, lamomega, theJPCLS, fitPair, fitParmOmMass, fitParmOmWidth, theData,1)+calcOmegaProdPartAmp(Minit, lamomega, theJPCLS, fitPair, fitParmOmMass, fitParmOmWidth, theData,2);
  }
    
  return result;

}

complex<double>  OmegaPiLhGammaBw::calcOmegaProdPartAmp(Spin Minit, Spin lamomega, boost::shared_ptr<const JPCLS> theJPCLS, pair<double, double> fitVal,  double fitParmOmMass, double fitParmOmWidth, OmegaPiData::OmPiEvtData* theData, unsigned int piGamCombi){

  complex<double> result(0.,0.);
  if (fabs(lamomega) > theJPCLS->J) return result;
  
  double theMag=fitVal.first;
  double thePhi=fitVal.second;
  complex<double> expiphi(cos(thePhi), sin(thePhi));
  
  

  if (piGamCombi==1){
    Vector4<double> omega4Vec(theData->omegaHeliCm4Vec.E(), theData->omegaHeliCm4Vec.Px(), theData->omegaHeliCm4Vec.Py(), theData->omegaHeliCm4Vec.Pz());
    result=sqrt(2*theJPCLS->L+1)*sqrt(3)*theMag*expiphi*Clebsch(theJPCLS->L,0,1, lamomega,theJPCLS->J, lamomega)*conj(theData->Dfp[theJPCLS->J][Minit][lamomega])*BreitWigner(omega4Vec, fitParmOmMass, fitParmOmWidth);
  }
  
  else if (piGamCombi==2){ 
    Vector4<double> omega4Vec=Vector4<double>(theData->omegaHeliCm4Vec2.E(), theData->omegaHeliCm4Vec2.Px(), theData->omegaHeliCm4Vec2.Py(), theData->omegaHeliCm4Vec2.Pz());
    result=sqrt(2*theJPCLS->L+1)*sqrt(3)*theMag*expiphi*Clebsch(theJPCLS->L,0,1, lamomega,theJPCLS->J, lamomega)*conj(theData->Dfp2[theJPCLS->J][Minit][lamomega])*BreitWigner(omega4Vec, fitParmOmMass, fitParmOmWidth);
    
  }
  
  else{
    Alert << "wrong combination in OmegaPiLhGammaBw::calcOmegaProdPartAmp" << endmsg;
    exit(1);
  }

  if (Minit==-1 && theJPCLS->P==-1 && theJPCLS->C==-1){
    int theJ=theJPCLS->J;
    if( theJ!=0 && (theJ%2==0 ) ){ //J=even => negativ sign, due to CGs of the pbar p production
      result=-result;
    }  
  }  

  return result;
}

complex<double> OmegaPiLhGammaBw::spinDensity(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal){
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmSinglet=theParamVal.omegaProdSinglet;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet0=theParamVal.omegaProdTriplet0;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet1=theParamVal.omegaProdTriplet1;

  double omegaMassFit=theParamVal.omegaMass.first;
  double omegaWidthFit=theParamVal.omegaMass.second;

  complex<double> MsingletAmpGM1=calcOmegaProdAmp(0, M, fitParmSinglet, omegaMassFit, omegaWidthFit, theData);  
  complex<double> Mtriplet0AmpGM1=calcOmegaProdAmp(0, M, fitParmTriplet0, omegaMassFit, omegaWidthFit, theData);
  complex<double> MtripletP1AmpGM1=calcOmegaProdAmp(1, M, fitParmTriplet1, omegaMassFit, omegaWidthFit, theData);
  complex<double> MtripletM1AmpGM1=calcOmegaProdAmp(-1, M, fitParmTriplet1, omegaMassFit, omegaWidthFit, theData);

  complex<double> M_singletAmpGM1=calcOmegaProdAmp(0, M_, fitParmSinglet, omegaMassFit, omegaWidthFit, theData);  
  complex<double> M_triplet0AmpGM1=calcOmegaProdAmp(0, M_, fitParmTriplet0, omegaMassFit, omegaWidthFit, theData);
  complex<double> M_tripletP1AmpGM1=calcOmegaProdAmp(1, M_, fitParmTriplet1, omegaMassFit, omegaWidthFit, theData);
  complex<double> M_tripletM1AmpGM1=calcOmegaProdAmp(-1, M_, fitParmTriplet1, omegaMassFit, omegaWidthFit, theData);

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
      MsingletAmpGM1=calcOmegaProdAmp(0,M1, fitParmSinglet, omegaMassFit, omegaWidthFit, theData);  
      Mtriplet0AmpGM1=calcOmegaProdAmp(0,M1, fitParmTriplet0, omegaMassFit, omegaWidthFit, theData);
      MtripletP1AmpGM1=calcOmegaProdAmp(1,M1, fitParmTriplet1, omegaMassFit, omegaWidthFit, theData);
      MtripletM1AmpGM1=calcOmegaProdAmp(-1,M1, fitParmTriplet1, omegaMassFit, omegaWidthFit, theData);
      theNorm += (norm(MsingletAmpGM1)+norm(Mtriplet0AmpGM1)+norm(MtripletP1AmpGM1)+norm(MtripletM1AmpGM1));
    }
  }


  return (result/theNorm);

}




void OmegaPiLhGammaBw::getFitParamVal(OmegaPiData::fitParamVal& theParamVal, const std::vector<double>& par) const{

  OmegaPiLhGamma::getFitParamVal(theParamVal, par);

  std::pair<double, double> omMass;

  omMass.first=par[par.size()-2];
  omMass.second=par[par.size()-1];

  theParamVal.omegaMass=omMass;
}

void OmegaPiLhGammaBw::setGenevaFitParamVal( boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObjectCollection> theGbdc_ptr ){

  OmegaPiLhGamma::setGenevaFitParamVal( theGbdc_ptr );

  boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObject> gbdMass_ptr(new Gem::Geneva::GConstrainedDoubleObject(0.7, 0.85) ); //omega mass
  theGbdc_ptr->push_back(gbdMass_ptr);

  boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObject> gbdWidth_ptr(new Gem::Geneva::GConstrainedDoubleObject(0.005, 0.04) ); //omega mass
  theGbdc_ptr->push_back(gbdWidth_ptr);

}


void OmegaPiLhGammaBw::setMnUsrParams(MnUserParameters& upar){

  OmegaPiLhGamma::setMnUsrParams(upar);

  upar.Add("omegaMass", 0.782, .001, 0.7, 0.85);
  upar.Add("omegaWidth", 0.008, .0005, 0.005, 0.04);
}

void OmegaPiLhGammaBw::setMnUsrParams(MnUserParameters& upar, OmegaPiData::fitParamVal &finalFitParm){
  
  OmegaPiLhGamma::setMnUsrParams(upar, finalFitParm);

  double omegaMassFit=finalFitParm.omegaMass.first;
  double omegaWidthFit=finalFitParm.omegaMass.second;
  
  upar.Add("omegaMass", omegaMassFit, .001, 0.7, 0.85);
  upar.Add("omegaWidth", omegaWidthFit, .0005, 0.005, 0.04);

}



void  OmegaPiLhGammaBw::setMnUsrParams(MnUserParameters& upar, minuitStartParam &theStartParam)
{
  OmegaPiLhGamma::setMnUsrParams(upar, theStartParam);

  std::string strName = "Omega";

  ParameterMap::iterator it;

  it = theStartParam.getParamMap().find(strName);
  if (it != theStartParam.getParamMap().end()){
    upar.Add("omegaMass", it->second[0], .001, 0.7, 0.85);
    upar.Add("omegaWidth", it->second[1], .0005, 0.005, 0.04);
    DebugMsg << "\nUsing user start Parameter for omega " << "with mass= " << it->second[0]
	   << " and width= " <<  it->second[1] << "\n"; 
  }
}

void OmegaPiLhGammaBw::dumpCurrentResult(std::ostream& os, const OmegaPiData::fitParamVal& fitParmVal) const{

  OmegaPiLhGamma::dumpCurrentResult(os, fitParmVal);

  pair<double, double>  theOmegaMass=fitParmVal.omegaMass;
  
  os << "Omega" << "\t" << theOmegaMass.first  << "\t" << theOmegaMass.second  << "\n"; 

}


void OmegaPiLhGammaBw::printFitParams(std::ostream& os, const OmegaPiData::fitParamVal& fitParmVal){

  OmegaPiLhGamma::printFitParams(os, fitParmVal);

  os << "***fit parameter for omega state*** " << "\n"; 
  os <<"\t mass:" << fitParmVal.omegaMass.first <<"\t width:" << fitParmVal.omegaMass.second  << "\n"; 
  
}

void OmegaPiLhGammaBw::print(std::ostream& os) const{
  os << "OmegaPiLhGammaBw\n";
}
