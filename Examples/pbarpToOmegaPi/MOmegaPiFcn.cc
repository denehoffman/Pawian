//#include <getopt.h>
//#include <fstream>
//#include <string>
#include <math.h>
#include <stdio.h>

#include "Minuit2/MnUserParameters.h"

#include "Examples/pbarpToOmegaPi/MOmegaPiFcn.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiLh.hh"
#include "ErrLogger/ErrLineLog.hh"

using namespace ROOT::Minuit2;

MOmegaPiFcn::MOmegaPiFcn(boost::shared_ptr<OmegaPiLh> omegaPiLh) :
  _omegaPiLhPtr(omegaPiLh),
  _barpToOmegaPi0States(omegaPiLh->omegaPi0States())
{
  if (0==_omegaPiLhPtr) ErrMsg(fatal) << "OmegaPiLh pointer is 0 !!!!" << endmsg; 
  
}

MOmegaPiFcn::~MOmegaPiFcn()
{
}

double MOmegaPiFcn::operator()(const std::vector<double>& par) const
{
  OmegaPiData::fitParamVal theFitParmValTmp;
  setFitParamVal(theFitParmValTmp, par);
 
  double result=_omegaPiLhPtr->calcLogLh(theFitParmValTmp);

  //  print fit paramss
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaSinglet=_barpToOmegaPi0States->jpclsSinglet();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0=_barpToOmegaPi0States->jpclsTriplet0();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1=_barpToOmegaPi0States->jpclsTriplet1();
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  ErrMsg(debugging) << "logLh= " << result <<endmsg;   
  ErrMsg(debugging) << "***fit parameter singlet states*** " <<endmsg;  
  for ( itJPCLS=JPCLSOmegaSinglet.begin(); itJPCLS!=JPCLSOmegaSinglet.end(); ++itJPCLS){
    ErrMsg(debugging)<< (*itJPCLS)->name()<< endmsg;
    std::pair<double, double> tmpParam=theFitParmValTmp.omegaProdSinglet[(*itJPCLS)];
    ErrMsg(debugging) <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }
  ErrMsg(debugging) << "***fit parameter triplet m=0 states*** " <<endmsg;  
  for ( itJPCLS=JPCLSOmegaTriplet0.begin(); itJPCLS!=JPCLSOmegaTriplet0.end(); ++itJPCLS){
    ErrMsg(debugging)<< (*itJPCLS)->name()<< endmsg;
    std::pair<double, double> tmpParam=theFitParmValTmp.omegaProdTriplet0[(*itJPCLS)];
    ErrMsg(debugging) <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }
  ErrMsg(debugging) << "***fit parameter triplet m=1 states*** " <<endmsg;  
  for ( itJPCLS=JPCLSOmegaTriplet1.begin(); itJPCLS!=JPCLSOmegaTriplet1.end(); ++itJPCLS){
    ErrMsg(debugging)<< (*itJPCLS)->name()<< endmsg;
    std::pair<double, double> tmpParam=theFitParmValTmp.omegaProdTriplet1[(*itJPCLS)];
    ErrMsg(debugging) <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }
    ErrMsg(debugging) << endmsg;  
  return result;
}

double MOmegaPiFcn::Up() const 
{
return .5;
}

void  MOmegaPiFcn::setMnUsrParams(MnUserParameters& upar){
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaSinglet=_barpToOmegaPi0States->jpclsSinglet();

  int counter=0;
  for ( itJPCLS=JPCLSOmegaSinglet.begin(); itJPCLS!=JPCLSOmegaSinglet.end(); ++itJPCLS){
    //now fill the fitParameterMap
    std::string magStr=(*itJPCLS)->name()+"S"+"mag";
    std::string phiStr=(*itJPCLS)->name()+"S"+"phi";
    upar.Add(magStr, 0.5, .1, -10., 10.);
    if (counter>0) upar.Add(phiStr, 0., .1, -M_PI, M_PI);
    counter++; 
  }

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0=_barpToOmegaPi0States->jpclsTriplet0();
  counter=0;
  for ( itJPCLS=JPCLSOmegaTriplet0.begin(); itJPCLS!=JPCLSOmegaTriplet0.end(); ++itJPCLS){
    //now fill the fitParameterMap
    std::string magStr=(*itJPCLS)->name()+"T0"+"mag";
    std::string phiStr=(*itJPCLS)->name()+"T0"+"phi";
    upar.Add(magStr, 0.5, .1, -10., 10.);
    if (counter>0) upar.Add(phiStr, 0., .1, -M_PI, M_PI);
    counter++; 
  }

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1=_barpToOmegaPi0States->jpclsTriplet1();
  counter=0;
  for ( itJPCLS=JPCLSOmegaTriplet1.begin(); itJPCLS!=JPCLSOmegaTriplet1.end(); ++itJPCLS){
    //now fill the fitParameterMap
    std::string magStr=(*itJPCLS)->name()+"T1"+"mag";
    std::string phiStr=(*itJPCLS)->name()+"T1"+"phi";
    upar.Add(magStr, 0.5, .1, -10., 10.);
    if (counter>0) upar.Add(phiStr, 0., .1, -M_PI, M_PI);
    counter++; 
  }


}


void MOmegaPiFcn::setFitParamVal(OmegaPiData::fitParamVal& theParamVal, const std::vector<double>& par) const{

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaSinglet=_barpToOmegaPi0States->jpclsSinglet();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0=_barpToOmegaPi0States->jpclsTriplet0();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1=_barpToOmegaPi0States->jpclsTriplet1();

  if (par.size()!= JPCLSOmegaSinglet.size()*2+JPCLSOmegaTriplet0.size()*2+JPCLSOmegaTriplet1.size()*2-3) ErrMsg(fatal) << "size of parameters wrong!!! par.size()=" << par.size() << "\tJPCLSOmegaSinglet.size()+JPCLSOmegaTriplet0.size()+JPCLSOmegaTriplet1.size()-3=" << JPCLSOmegaSinglet.size()*2+JPCLSOmegaTriplet0.size()*2+JPCLSOmegaTriplet1.size()*2-3 << endmsg;

  int counter=0;
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

