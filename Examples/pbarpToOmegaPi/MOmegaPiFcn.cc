//#include <getopt.h>
//#include <fstream>
//#include <string>
#include <math.h>
#include <stdio.h>
#include "Minuit2/MnPrint.h"

#include "Minuit2/MnUserParameters.h"

#include "Examples/pbarpToOmegaPi/MOmegaPiFcn.hh"
#include "Examples/pbarpToOmegaPi/AbsOmegaPiLh.hh"
#include "ErrLogger/ErrLogger.hh"

using namespace ROOT::Minuit2;

MOmegaPiFcn::MOmegaPiFcn(boost::shared_ptr<AbsOmegaPiLh> absOmegaPiLh) :
  _omegaPiLhPtr(absOmegaPiLh),
  _barpToOmegaPi0States(absOmegaPiLh->omegaPi0States())
{
  if (0==_omegaPiLhPtr) { Alert << "OmegaPiLh pointer is 0 !!!!" << endmsg; exit(1); }
  
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

/*  DebugMsg << "logLh= " << result <<endmsg;   
  DebugMsg << "***fit parameter singlet states*** " <<endmsg;  
  for ( itJPCLS=JPCLSOmegaSinglet.begin(); itJPCLS!=JPCLSOmegaSinglet.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< endmsg;
    std::pair<double, double> tmpParam=theFitParmValTmp.omegaProdSinglet[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }
  DebugMsg << "***fit parameter triplet m=0 states*** " <<endmsg;  
  for ( itJPCLS=JPCLSOmegaTriplet0.begin(); itJPCLS!=JPCLSOmegaTriplet0.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< endmsg;
    std::pair<double, double> tmpParam=theFitParmValTmp.omegaProdTriplet0[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }
  DebugMsg << "***fit parameter triplet m=1 states*** " <<endmsg;  
  for ( itJPCLS=JPCLSOmegaTriplet1.begin(); itJPCLS!=JPCLSOmegaTriplet1.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< endmsg;
    std::pair<double, double> tmpParam=theFitParmValTmp.omegaProdTriplet1[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }
  DebugMsg << endmsg;*/  
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

    upar.Add(magStr, 0.5, .1, 0., 2.);
    if (counter>0) upar.Add(phiStr, 0., .1, -3.*M_PI, 3.*M_PI);
    counter++; 
  }

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0=_barpToOmegaPi0States->jpclsTriplet0();
  counter=0;
  for ( itJPCLS=JPCLSOmegaTriplet0.begin(); itJPCLS!=JPCLSOmegaTriplet0.end(); ++itJPCLS){
    //now fill the fitParameterMap
    std::string magStr=(*itJPCLS)->name()+"T0"+"mag";
    std::string phiStr=(*itJPCLS)->name()+"T0"+"phi";

    upar.Add(magStr, 0.5, .1, 0., 2.);
    if (counter>0) upar.Add(phiStr, 0., .1, -3.*M_PI, 3.*M_PI);
    counter++; 
  }

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1=_barpToOmegaPi0States->jpclsTriplet1();
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

inline void FillUpars(MnUserParameters& upar, 
                      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmS,
                      const string &theSuffix
                     )
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

      upar.Add(magStr, theMag, .1, 0., 1.);
      if (!bFirst) upar.Add(phiStr, thePhi, .1, -3.*M_PI, 3.*M_PI);
      if(bFirst) bFirst=false;
   } 
}
 
void  MOmegaPiFcn::setMnUsrParams(MnUserParameters& upar, OmegaPiData::fitParamVal &finalFitParm)
{

   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmSinglet=finalFitParm.omegaProdSinglet;
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet0=finalFitParm.omegaProdTriplet0;
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > fitParmTriplet1=finalFitParm.omegaProdTriplet1;
   FillUpars(upar,fitParmSinglet,"S");
   FillUpars(upar,fitParmTriplet0,"T0");
   FillUpars(upar,fitParmTriplet1,"T1");
}

//This Method fills the Minuit User Parameters with parameters for given state and if supllied initializes
//the parameter with start parameters given by user in theStartParam map.
inline void FillUpars(MnUserParameters& upar, 
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

void  MOmegaPiFcn::setMnUsrParams(MnUserParameters& upar, minuitStartParam &theStartParam)
{
  FillUpars(upar,theStartParam,_barpToOmegaPi0States->jpclsSinglet(),"S");
  FillUpars(upar,theStartParam,_barpToOmegaPi0States->jpclsTriplet0(),"T0");
  FillUpars(upar,theStartParam,_barpToOmegaPi0States->jpclsTriplet1(),"T1");
}

void MOmegaPiFcn::setFitParamVal(OmegaPiData::fitParamVal& theParamVal, const std::vector<double>& par) const{

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaSinglet=_barpToOmegaPi0States->jpclsSinglet();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet0=_barpToOmegaPi0States->jpclsTriplet0();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSOmegaTriplet1=_barpToOmegaPi0States->jpclsTriplet1();

  if (par.size()!= JPCLSOmegaSinglet.size()*2+JPCLSOmegaTriplet0.size()*2+JPCLSOmegaTriplet1.size()*2-3) {
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

