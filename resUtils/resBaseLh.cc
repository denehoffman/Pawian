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

// resBaseLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "resUtils/resBaseLh.hh"
#include "resUtils/resEnv.hh"
#include "resUtils/resReaction.hh"
#include "resUtils/resEnv.hh"
#include "PwaUtils/HeliDecAmps.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


resBaseLh::resBaseLh() :
  AbsLh(resEnv::instance())
  ,_highestJFsp(0)
  ,_isHighestJaPhoton(true)
{
  initialize();
}


resBaseLh::~resBaseLh()
{;
}


complex<double> resBaseLh::calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData){

  complex<double> result(0.,0.);

  return result;
}




complex<double> resBaseLh::calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData, 
					     std::map <std::shared_ptr<const JPCLS>,
					     std::vector< std::shared_ptr<AbsXdecAmp> >,
					     pawian::Collection::SharedPtrLess > pbarpAmps){
   complex<double> resultAmp(0.,0.);

   return resultAmp;
}


double resBaseLh::calcEvtIntensity(EvtData* theData, fitParams& theParamVal){

  double result=0.;

  int lamSteps=1;
  if(_isHighestJaPhoton) lamSteps=2;

  for (int lamHigestJFsp=-_highestJFsp; lamHigestJFsp<=_highestJFsp; lamHigestJFsp=lamHigestJFsp+lamSteps){

    for( Spin itLam=-_Jmother; itLam<=_Jmother; ++itLam){
      complex<double> lamItAmp(0.,0.);
      
      std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
      for( itDec=_decAmps.begin(); itDec!=_decAmps.end(); ++itDec){
	complex<double> currentDecAmp=(*itDec)->XdecAmp(itLam, theData, lamHigestJFsp);
	lamItAmp+=currentDecAmp;
      }
      result += norm(lamItAmp);
    }
  }

  if(_usePhasespace) result+=theParamVal.otherParams[_phasespaceKey];
  return result;  
  
}

void resBaseLh::print(std::ostream& os) const{

}


void  resBaseLh::initialize(){

  _Jmother = Spin(resEnv::instance()->motherParticle()->twoJ(), 2);

  std::vector<Particle*> fsParticles=resEnv::instance()->finalStateParticles();
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
  
  _resReactionPtr = resEnv::instance()->reaction();

  std::vector< std::shared_ptr<IsobarHeliDecay> > theDecs = _resReactionPtr->productionHeliDecays();
  std::vector< std::shared_ptr<IsobarHeliDecay> >::iterator it;
  for (it=theDecs.begin(); it!=theDecs.end(); ++it){
    //    std::shared_ptr<AbsDecay> currentDec((*it).get() );
    //    std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(currentDec);
    std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp((*it)->absDecPtr());
    _decAmps.push_back(currentAmp);
  }

}




