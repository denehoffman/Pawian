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

// pbarpTensorLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "pbarpUtils/pbarpTensorLh.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "pbarpUtils/pbarpReaction.hh"
#include "PwaUtils/LSDecAmps.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


pbarpTensorLh::pbarpTensorLh() :
  pbarpBaseLh()
{
  initialize();
}


pbarpTensorLh::~pbarpTensorLh()
{;
}


void pbarpTensorLh::print(std::ostream& os) const{

}


void  pbarpTensorLh::initialize(){

  std::vector< std::shared_ptr<IsobarTensorDecay> > theDecs = _pbarpReactionPtr->productionTensorDecays();
  std::vector< std::shared_ptr<IsobarTensorDecay> >::iterator it;
  for (it=theDecs.begin(); it!=theDecs.end(); ++it){
    std::shared_ptr<AbsDecay> currentDec( (*it).get() );
    std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(currentDec);
    _decAmps.push_back(currentAmp);
  }

  std::vector< std::shared_ptr<const JPCLS> > jpclsSingletStates=_pbarpReactionPtr->jpclsSingletStates();
  fillMap(jpclsSingletStates, _decAmps, _decAmpsSinglet);

  std::vector< std::shared_ptr<const JPCLS> > jpclsTriplet0States=_pbarpReactionPtr->jpclsTriplet0States();
  fillMap(jpclsTriplet0States, _decAmps, _decAmpsTriplet0);

  std::vector< std::shared_ptr<const JPCLS> > jpclsTripletp1States=_pbarpReactionPtr->jpclsTripletp1States();
  fillMap(jpclsTripletp1States, _decAmps, _decAmpsTripletp1);

  std::vector< std::shared_ptr<const JPCLS> > jpclsTripletm1States=_pbarpReactionPtr->jpclsTripletm1States();
  fillMap(jpclsTripletm1States, _decAmps, _decAmpsTripletm1); 
  fillIsos();

  double pbarMass = _absEnv->particleTable()->particle("antiproton")->mass();
  double pMass = _absEnv->particleTable()->particle("proton")->mass();
  double pbarMom = pbarpEnv::instance()->pbarMomentum();

  Vector4<double> pbar4Vec(sqrt(pbarMass*pbarMass+pbarMom*pbarMom), 0.,0., pbarMom);
  Vector4<double> p4Vec(pMass, 0.,0.,0.);
  Vector4<double> allVec=pbar4Vec+p4Vec;

  std::vector< std::shared_ptr<const JPCLS> > jpclsStatesAll=_pbarpReactionPtr->jpclsStates();
  std::vector< std::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  for(itJPCLS=jpclsStatesAll.begin(); itJPCLS!=jpclsStatesAll.end(); ++itJPCLS){
    std::shared_ptr<OrbitalTensor> currentTensorPtr(new OrbitalTensor((*itJPCLS)->L));
    currentTensorPtr->SetP4(pbar4Vec,p4Vec); 
    _orbTensorMap[(*itJPCLS)->L]=currentTensorPtr;

    std::shared_ptr<PolVector> currentPolVecPtr(new PolVector((*itJPCLS)->J));
    currentPolVecPtr->SetP4(allVec, allVec.M());
    _polVectorMap[(*itJPCLS)->J]=currentPolVecPtr; 
    }
}


