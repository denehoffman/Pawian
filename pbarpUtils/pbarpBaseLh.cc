//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

// pbarpBaseLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "pbarpUtils/pbarpBaseLh.hh"
#include "pbarpUtils/pbarpReaction.hh"
#include "pbarpUtils/PbarpChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/FsParticleProjections.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


pbarpBaseLh::pbarpBaseLh(ChannelID channelID) :
  AbsLh(channelID)
  ,_highestJFsp(0)
  ,_isHighestJaPhoton(true)
{
  initialize();
}


pbarpBaseLh::~pbarpBaseLh()
{;
}


complex<double> pbarpBaseLh::calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData, int J){

   complex<double> result(0.,0.);

   complex<double> singletAmpM1   = calcProdPartAmp( 0, M1, nameDec, theData, _decAmpsSinglet);
   complex<double> triplet0AmpM1  = calcProdPartAmp( 0, M1, nameDec, theData, _decAmpsTriplet0);
   complex<double> triplet1AmpM1  = calcProdPartAmp( 1, M1, nameDec, theData, _decAmpsTripletp1);
   complex<double> tripletm1AmpM1 = calcProdPartAmp(-1, M1, nameDec, theData, _decAmpsTripletm1);

   complex<double> singletAmpM2   = calcProdPartAmp( 0, M2, nameDec, theData, _decAmpsSinglet);
   complex<double> triplet0AmpM2  = calcProdPartAmp( 0, M2, nameDec, theData, _decAmpsTriplet0);
   complex<double> triplet1AmpM2  = calcProdPartAmp( 1, M2, nameDec, theData, _decAmpsTripletp1);
   complex<double> tripletm1AmpM2 = calcProdPartAmp(-1, M2, nameDec, theData, _decAmpsTripletm1);

   result = singletAmpM1   * conj(singletAmpM2) +
            triplet0AmpM1  * conj(triplet0AmpM2) +
            triplet1AmpM1  * conj(triplet1AmpM2) +
            tripletm1AmpM1 * conj(tripletm1AmpM2);

   double norm = calcSpinDensityNorm(nameDec, theData, J);

   return (result / norm);
}



double pbarpBaseLh::calcSpinDensityNorm(std::string& nameDec, EvtData* theData, int J){

   double result = 0;

   for (Spin M=-J; M<=J; ++M){
     complex<double> singletAmpM1   = calcProdPartAmp( 0, M, nameDec, theData, _decAmpsSinglet);
     complex<double> triplet0AmpM1  = calcProdPartAmp( 0, M, nameDec, theData, _decAmpsTriplet0);
     complex<double> triplet1AmpM1  = calcProdPartAmp( 1, M, nameDec, theData, _decAmpsTripletp1);
     complex<double> tripletm1AmpM1 = calcProdPartAmp(-1, M, nameDec, theData, _decAmpsTripletm1);

     result += (norm(singletAmpM1) + norm(triplet0AmpM1) + norm(triplet1AmpM1) + norm(tripletm1AmpM1));
   }

   return result;
}


complex<double> pbarpBaseLh::calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData,
					     std::map <std::shared_ptr<const JPCLS>,
					     std::vector< std::shared_ptr<AbsXdecAmp> >,
					     pawian::Collection::SharedPtrLess > pbarpAmps){

  std::vector< std::shared_ptr<AbsXdecAmp> >::iterator itDecAll;
  for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
    (*itDecAll)->calcDynamics(theData);
  }

   complex<double> resultAmp(0.,0.);

   std::map <std::shared_ptr<const JPCLS>,
             std::vector< std::shared_ptr<AbsXdecAmp> >,
             pawian::Collection::SharedPtrLess >::iterator it;

   for(it=pbarpAmps.begin(); it!=pbarpAmps.end(); ++it){
     complex<double> tmpAmp(0.,0.);
     std::shared_ptr<const JPCLS> theJPCLS=it->first;
     double iso1Val=_currentParamJPCIsos1[theJPCLS];
     double iso0Val=_currentParamJPCIsos0[theJPCLS];

     std::vector<std::shared_ptr<AbsXdecAmp> > decAmps=it->second;
     std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       Particle* particle1 = (*itDec)->absDec()->daughter1Part();
       Particle* particle2 = (*itDec)->absDec()->daughter2Part();

       double isoFactor=0;
       if((*itDec)->absDec()->motherIGJPC()->I==1) isoFactor=iso1Val;
       else isoFactor=iso0Val;

       Spin lamFs=0;
       if(particle1->name() == nameDec)
	 tmpAmp+=isoFactor * (*itDec)->XdecPartAmp(lamX, lamDec, 1,  theData, lamFs);
       else if(particle2->name() == nameDec)
	 tmpAmp+=isoFactor * (*itDec)->XdecPartAmp(lamX, lamDec, 2,  theData, lamFs);
     }

     double theMag=_currentParamMags[theJPCLS];
     double thePhi=_currentParamPhis[theJPCLS];
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=theJPCLS->preFactor*theMag*expi;
     resultAmp+=tmpAmp;
   }

   return resultAmp;
}


double pbarpBaseLh::calcEvtIntensity(EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar){

  double result=0.;

  std::vector< std::shared_ptr<AbsXdecAmp> >::iterator itDecAll;
  for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
    (*itDecAll)->calcDynamics(theData);
  }

  std::map <std::shared_ptr<const JPCLS>, std::vector< std::shared_ptr<AbsXdecAmp> >, pawian::Collection::SharedPtrLess >::iterator it;

  std::vector< std::vector<Spin> > spinProjections=_fsParticleProjections->spinProjections();

  for (unsigned int projId=0; projId<spinProjections.size(); ++projId){
    for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
      (*itDecAll)->setSpinProjections(projId);
    }
    
    complex<double> singletAmp(0.,0.);
    Spin lampbarp=0;
    for(it=_decAmpsSinglet.begin(); it!=_decAmpsSinglet.end(); ++it){
      complex<double> tmpAmp(0.,0.);
      const std::shared_ptr<const JPCLS>& theJPCLS=it->first;
      double iso1Val=_currentParamJPCIsos1.at(theJPCLS);
      double iso0Val=_currentParamJPCIsos0.at(theJPCLS);

      std::vector<std::shared_ptr<AbsXdecAmp> >& decAmps=it->second;
      std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
      for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
	
	complex<double> currentDecAmp=(*itDec)->XdecAmp(lampbarp, theData);

	double isoFactor=0;
	if((*itDec)->absDec()->motherIGJPC()->I==1) isoFactor=iso1Val;
	else isoFactor=iso0Val;
	
	tmpAmp+= isoFactor*currentDecAmp;
      }
      
      tmpAmp*=theJPCLS->preFactor*_currentParamMagExpi.at(theJPCLS);
      singletAmp+=tmpAmp;
    }

    complex<double> triplet0Amp(0.,0.);
    lampbarp=0;
    for(it=_decAmpsTriplet0.begin(); it!=_decAmpsTriplet0.end(); ++it){
      complex<double> tmpAmp(0.,0.);
      const std::shared_ptr<const JPCLS>& theJPCLS=it->first;
      double iso1Val=_currentParamJPCIsos1.at(theJPCLS);
      double iso0Val=_currentParamJPCIsos0.at(theJPCLS);

      std::vector<std::shared_ptr<AbsXdecAmp> >& decAmps=it->second;
      std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
      for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
	complex<double> currentDecAmp=(*itDec)->XdecAmp(lampbarp, theData);

	double isoFactor=0;
	if((*itDec)->absDec()->motherIGJPC()->I==1) isoFactor=iso1Val;
	else isoFactor=iso0Val;
	
	tmpAmp+=isoFactor*currentDecAmp;
      }
      
      tmpAmp*=theJPCLS->preFactor*_currentParamMagExpi.at(theJPCLS);
      triplet0Amp+=tmpAmp;
    }

    complex<double> tripletp1Amp(0.,0.);
    lampbarp=1;
    for(it=_decAmpsTripletp1.begin(); it!=_decAmpsTripletp1.end(); ++it){
      complex<double> tmpAmp(0.,0.);
      const std::shared_ptr<const JPCLS>& theJPCLS=it->first;
      double iso1Val=_currentParamJPCIsos1.at(theJPCLS);
      double iso0Val=_currentParamJPCIsos0.at(theJPCLS);

      std::vector<std::shared_ptr<AbsXdecAmp> >& decAmps=it->second;
      std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
      for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
	complex<double> currentDecAmp=(*itDec)->XdecAmp(lampbarp, theData);

	double isoFactor=0;
	if((*itDec)->absDec()->motherIGJPC()->I==1) isoFactor=iso1Val;
	else isoFactor=iso0Val;
	
	tmpAmp+=isoFactor*currentDecAmp;
      }
      
      tmpAmp*=theJPCLS->preFactor*_currentParamMagExpi.at(theJPCLS);
      tripletp1Amp+=tmpAmp;
    }
    
    complex<double> tripletm1Amp(0.,0.);
    lampbarp=-1;
    for(it=_decAmpsTripletm1.begin(); it!=_decAmpsTripletm1.end(); ++it){
      complex<double> tmpAmp(0.,0.);
      const std::shared_ptr<const JPCLS>& theJPCLS=it->first;
      double iso1Val=_currentParamJPCIsos1.at(theJPCLS);
      double iso0Val=_currentParamJPCIsos0.at(theJPCLS);

      std::vector<std::shared_ptr<AbsXdecAmp> >& decAmps=it->second;
      std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
      for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
  	complex<double> currentDecAmp=(*itDec)->XdecAmp(lampbarp, theData);

        double isoFactor=0;
  	if((*itDec)->absDec()->motherIGJPC()->I==1) isoFactor=iso1Val;
  	else isoFactor=iso0Val;
	
  	tmpAmp+=isoFactor*currentDecAmp;
      }
      
      tmpAmp*=theJPCLS->preFactor*_currentParamMagExpi.at(theJPCLS);
      tripletm1Amp+=tmpAmp;
    }

    result += norm(singletAmp)+ norm(triplet0Amp)+ norm(tripletp1Amp)+ norm(tripletm1Amp);
  }
  
  if(_usePhasespace) result+=fitPar->Value(_phasespaceKey);
  
  result *= fitPar->Value(_channelScaleParam);
  
  return result;
  
}

void pbarpBaseLh::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){

  AbsLh::fillDefaultParams(fitPar);

  std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess>::iterator itIso;
  for(itIso=_currentParamJPCIsos1.begin(); itIso!=_currentParamJPCIsos1.end(); ++itIso){
    fitPar->Add("Iso1"+(*itIso).first->name(), (*itIso).second, 0.5);
    fitPar->SetLimits("Iso1"+(*itIso).first->name(), 0., 1.);
   }

  double magFactor=1./sqrt(_jpclsStates.size());
  std::vector< std::shared_ptr<const JPCLS> >::iterator it;
  for ( it = _jpclsStates.begin(); it!=_jpclsStates.end(); ++it){
    //fill magnitude
    std::string magName=(*it)->name()+"pbarp"+"Mag";
    double valMag=magFactor;
    double errMag=magFactor/2.;
    //    double minMag=magFactor-6.*errMag;
    // if (minMag<0.) minMag=0.;
    // double maxMag=magFactor+30.*errMag;
    fitPar->Add(magName, valMag, errMag);
    //    fitPar->SetLimits(magName, minMag, maxMag);

    std::string phiName=(*it)->name()+"pbarp"+"Phi";
    double valPhi=0.;
    double errPhi=0.2;
    //no limits for phi parameter
    fitPar->Add(phiName, valPhi, errPhi);
  }
  
}

void pbarpBaseLh::print(std::ostream& os) const{

}


void  pbarpBaseLh::initialize(){

  std::vector<Particle*> fsParticles=std::static_pointer_cast<PbarpChannelEnv>(GlobalEnv::instance()->PbarpChannel(_channelID))->finalStateParticles();
  std::vector<Particle*>::iterator itParticle;
  bool highestJFound=false;

  for (itParticle=fsParticles.begin(); itParticle != fsParticles.end(); ++itParticle){
    int current2J = (*itParticle)->twoJ();
    if(current2J>0){
      if(highestJFound){
	Alert << "final states with more than 1 particles with J>0 not supported!!!!" << endmsg;
	exit(1);
      }
      _highestJFsp=int(current2J/2);
      if( (*itParticle)->name() != "photon" ) _isHighestJaPhoton=false;
    }
  }

  _pbarpReactionPtr = std::static_pointer_cast<PbarpChannelEnv>(GlobalEnv::instance()->PbarpChannel(_channelID))->reaction();
  _jpclsStates=_pbarpReactionPtr->jpclsStates();
  _igjpcStates=_pbarpReactionPtr->igjpcStates();
}

void pbarpBaseLh::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  AbsLh::updateFitParams(fitPar);

  std::map< std::shared_ptr<const jpcRes>, double, pawian::Collection::SharedPtrLess>::iterator itIso;
  for(itIso=_currentParamJPCIsos1.begin(); itIso!=_currentParamJPCIsos1.end(); ++itIso){
    double theVal=fitPar->Value("Iso1"+(*itIso).first->name());   
    (*itIso).second=theVal;
    _currentParamJPCIsos0[(*itIso).first]=sqrt(1.-theVal*theVal);
   }

  std::vector< std::shared_ptr<const JPCLS> >::iterator it;
  for ( it = _jpclsStates.begin(); it!=_jpclsStates.end(); ++it){
    std::string magName=(*it)->name()+"pbarp"+"Mag";
    std::string phiName=(*it)->name()+"pbarp"+"Phi";
    double theMag=fabs(fitPar->Value(magName));
    double thePhi=fitPar->Value(phiName);

    _currentParamMags[*it]=theMag;
    _currentParamPhis[*it]=thePhi; 

    complex<double> expi(cos(thePhi), sin(thePhi));
    _currentParamMagExpi[*it]=theMag*expi;
   }
}

void pbarpBaseLh::fillMap(std::vector< std::shared_ptr<const JPCLS> >& pbarpLSs, std::vector<std::shared_ptr<AbsXdecAmp> >& decs, std::map< std::shared_ptr<const JPCLS>, std::vector<std::shared_ptr<AbsXdecAmp> >, pawian::Collection::SharedPtrLess >& toFill){


  std::vector< std::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  for (itJPCLS = pbarpLSs.begin(); itJPCLS != pbarpLSs.end(); ++itJPCLS){
    std::vector<std::shared_ptr<AbsXdecAmp> > currentAmpVec;

    std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itAmp;
    for(itAmp=decs.begin(); itAmp!=decs.end(); ++itAmp){
      if( (*(*itAmp)->jpcPtr())==(*(*itJPCLS)) )  currentAmpVec.push_back(*itAmp);
    }
    toFill[(*itJPCLS)]=currentAmpVec;
  }

}

void pbarpBaseLh::fillIsos(){
  //first look for iso0 decay amplitudes and fir it in the map
  std::vector< std::shared_ptr<AbsXdecAmp> >::iterator it;
  for(it= _decAmps.begin(); it!= _decAmps.end(); ++it){
    if((*it)->absDec()->motherIGJPC()->I==0){
      _iso0DecAmps.push_back(*it);
    }
    if((*it)->absDec()->motherIGJPC()->I==1){
      _iso1DecAmps.push_back(*it);
    }
  }


  for(it= _iso0DecAmps.begin(); it!= _iso0DecAmps.end(); ++it){
    _currentParamJPCIsos0[(*it)->jpcPtr()]=1.;
    _currentParamJPCIsos1[(*it)->jpcPtr()]=0.;
  }

  for(it= _iso1DecAmps.begin(); it!= _iso1DecAmps.end(); ++it){
    _currentParamJPCIsos1[(*it)->jpcPtr()]=1.;
    _currentParamJPCIsos0[(*it)->jpcPtr()]=0.;
  }

  //now look if iso0 and iso1 exisist
  for(it= _iso0DecAmps.begin(); it!= _iso0DecAmps.end(); ++it){
    // find corresponding iso1 partner
    std::vector< std::shared_ptr<AbsXdecAmp> >::iterator itIso1;
    for(itIso1= _iso1DecAmps.begin(); itIso1!= _iso1DecAmps.end(); ++itIso1){
      if( (*(*it)->jpcPtr()) == (*(*itIso1)->jpcPtr()) ){
	_currentParamJPCIsos0[(*it)->jpcPtr()]=sqrt(.5);
	_currentParamJPCIsos1[(*it)->jpcPtr()]=sqrt(.5);
      }
    }

  }


  //now look if iso0 and iso1 exisist for the same amplitude
  for(it= _iso0DecAmps.begin(); it!= _iso0DecAmps.end(); ++it){
    std::vector< std::shared_ptr<AbsXdecAmp> >::iterator it2;
    for(it2= _iso1DecAmps.begin(); it2!= _iso1DecAmps.end(); ++it2){
      if((*it)->jpcDecsName()==(*it2)->jpcDecsName()){
	//found
	std::pair<std::string, std::string> currentNamePair=make_pair((*it)->name(), (*it2)->name());
	_iso0Iso1NameMap[(*it)->isoKey()]=currentNamePair;
	_currentParamJPCIsos1[(*it)->jpcPtr()]=sqrt(.7);
	_currentParamJPCIsos0[(*it)->jpcPtr()]=sqrt(.3);
      }
    }
  }
}
