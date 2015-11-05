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
#include "pbarpUtils/pbarpReaction.hh"
#include "pbarpUtils/PbarpChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "ErrLogger/ErrLogger.hh"
#include "FitParams/AbsPawianParameters.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


pbarpTensorLh::pbarpTensorLh(ChannelID channelID) :
   pbarpBaseLh(channelID)
{
  initialize();
}


pbarpTensorLh::~pbarpTensorLh()
{;
}


void pbarpTensorLh::print(std::ostream& os) const{

}

double pbarpTensorLh::calcEvtIntensity(EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar){

  double result=0.;

  std::vector< std::shared_ptr<AbsXdecAmp> >::iterator itDecAll;
  for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
    (*itDecAll)->calcDynamics(theData);
  }

  std::map <std::shared_ptr<const JPCLS>, std::vector< std::shared_ptr<AbsXdecAmp> >, pawian::Collection::SharedPtrLess >::iterator it;

  Spin lamSteps=1;
  if(_isHighestJaPhoton) lamSteps=2;

  for (Spin lamHigestJFsp=-_highestJFsp; lamHigestJFsp<=_highestJFsp; lamHigestJFsp=lamHigestJFsp+lamSteps){

  complex<double> singletAmp(0.,0.);
  Spin lampbarp=0;
  for(it=_decAmpsSinglet.begin(); it!=_decAmpsSinglet.end(); ++it){
    complex<double> tmpAmp(0.,0.);
    const std::shared_ptr<const JPCLS>& theJPCLS=it->first;
    double iso1Val=_currentParamJPCIsos1[theJPCLS];
    double iso0Val=_currentParamJPCIsos0[theJPCLS];

    std::vector<std::shared_ptr<AbsXdecAmp> >& decAmps=it->second;
     std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){

       complex<double> currentDecAmp=(*itDec)->XdecAmp(lampbarp, theData, lamHigestJFsp);

       double isoFactor=0;
       if((*itDec)->absDec()->motherIGJPC()->I==1) isoFactor=iso1Val;
       else isoFactor=iso0Val;

       tmpAmp+= isoFactor*currentDecAmp;
     }

     double theMag=_currentParamMags.at(theJPCLS);
     double thePhi=_currentParamPhis.at(theJPCLS);
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=_pbarpAmpSingletMap.at(theJPCLS)*theMag*expi;
     singletAmp+=tmpAmp;
  }

  complex<double> triplet0Amp(0.,0.);
  lampbarp=0;
  for(it=_decAmpsTriplet0.begin(); it!=_decAmpsTriplet0.end(); ++it){
    complex<double> tmpAmp(0.,0.);
    const std::shared_ptr<const JPCLS>& theJPCLS=it->first;
    double iso1Val=_currentParamJPCIsos1[theJPCLS];
    double iso0Val=_currentParamJPCIsos0[theJPCLS];

    std::vector<std::shared_ptr<AbsXdecAmp> >& decAmps=it->second;
     std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       complex<double> currentDecAmp=(*itDec)->XdecAmp(lampbarp, theData, lamHigestJFsp);

       double isoFactor=0;
       if((*itDec)->absDec()->motherIGJPC()->I==1) isoFactor=iso1Val;
       else isoFactor=iso0Val;

       tmpAmp+=isoFactor*currentDecAmp;
     }

     double theMag=_currentParamMags.at(theJPCLS);
     double thePhi=_currentParamPhis.at(theJPCLS);
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=_pbarpAmpTriplet0Map.at(theJPCLS)*theMag*expi;
     triplet0Amp+=tmpAmp;
  }

   complex<double> tripletp1Amp(0.,0.);
   lampbarp=1;
  for(it=_decAmpsTripletp1.begin(); it!=_decAmpsTripletp1.end(); ++it){
    complex<double> tmpAmp(0.,0.);
    const std::shared_ptr<const JPCLS>& theJPCLS=it->first;

    // std::cout << "\ntripletp1:" << std::endl;
    //  theJPCLS->print(std::cout);
    //  std::cout << "\t" << _pbarpAmpTripletp1Map.at(theJPCLS) << std::endl; 
   
     double iso1Val=_currentParamJPCIsos1[theJPCLS];
    double iso0Val=_currentParamJPCIsos0[theJPCLS];

    std::vector<std::shared_ptr<AbsXdecAmp> >& decAmps=it->second;
     std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
       complex<double> currentDecAmp=(*itDec)->XdecAmp(lampbarp, theData, lamHigestJFsp);

       double isoFactor=0;
       if((*itDec)->absDec()->motherIGJPC()->I==1) isoFactor=iso1Val;
       else isoFactor=iso0Val;

       tmpAmp+=isoFactor*currentDecAmp;
     }

     double theMag=_currentParamMags.at(theJPCLS);
     double thePhi=_currentParamPhis.at(theJPCLS);
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=_pbarpAmpTripletp1Map.at(theJPCLS)*theMag*expi;
     tripletp1Amp+=tmpAmp;
  }

   complex<double> tripletm1Amp(0.,0.);
   lampbarp=-1;
  for(it=_decAmpsTripletm1.begin(); it!=_decAmpsTripletm1.end(); ++it){
     complex<double> tmpAmp(0.,0.);
     const std::shared_ptr<const JPCLS> theJPCLS=it->first;
     //    std::cout << "\ntripletm1:" << std::endl;
     //     theJPCLS->print(std::cout);
     //     std::cout << "\t" << _pbarpAmpTripletm1Map.at(theJPCLS) << std::endl; 

     double iso1Val=_currentParamJPCIsos1[theJPCLS];
     double iso0Val=_currentParamJPCIsos0[theJPCLS];

     std::vector<std::shared_ptr<AbsXdecAmp> >& decAmps=it->second;
     std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
     for( itDec=decAmps.begin(); itDec!=decAmps.end(); ++itDec){
	complex<double> currentDecAmp=(*itDec)->XdecAmp(lampbarp, theData, lamHigestJFsp);

        double isoFactor=0;
	if((*itDec)->absDec()->motherIGJPC()->I==1) isoFactor=iso1Val;
	else isoFactor=iso0Val;

	tmpAmp+=isoFactor*currentDecAmp;
     }

     double theMag=_currentParamMags.at(theJPCLS);
     double thePhi=_currentParamPhis.at(theJPCLS);
     complex<double> expi(cos(thePhi), sin(thePhi));
     tmpAmp*=_pbarpAmpTripletm1Map.at(theJPCLS)*theMag*expi;
     tripletm1Amp+=tmpAmp;
  }

  result += 2.*norm(singletAmp)+ 2.*norm(triplet0Amp)+ norm(tripletp1Amp)+ norm(tripletm1Amp);
  }

  if(_usePhasespace) result+=fitPar->Value(_phasespaceKey);

  result *= fitPar->Value(_channelScaleParam);

  return result;

}


void  pbarpTensorLh::initialize(){

  std::vector< std::shared_ptr<IsobarTensorDecay> > theDecs = _pbarpReactionPtr->productionTensorDecays();
  std::vector< std::shared_ptr<IsobarTensorDecay> >::iterator it;
  for (it=theDecs.begin(); it!=theDecs.end(); ++it){
    std::shared_ptr<AbsDecay> currentDec( (*it).get() );
    std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(_channelID, currentDec);
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

  double pbarMass = GlobalEnv::instance()->particleTable()->particle("antiproton")->mass();
  double pMass = GlobalEnv::instance()->particleTable()->particle("proton")->mass();
  double pbarMom = std::static_pointer_cast<PbarpChannelEnv>(GlobalEnv::instance()->PbarpChannel(_channelID))->pbarMomentum();

  if(pbarMom>400.){
    pbarMass*=1000.;
    pMass*=1000.;
  }

  Vector4<double> pbar4Vec(sqrt(pbarMass*pbarMass+pbarMom*pbarMom), 0.,0., pbarMom);
  Vector4<double> p4Vec(pMass, 0.,0.,0.);
  Vector4<double> allVec=pbar4Vec+p4Vec;

  Tensor<complex<double> > cancelTerm=(pbar4Vec-p4Vec)/(pbarMass+pMass+allVec.M());

  DiracAntiSpinor pbarAntiSpinor; // pbar anti spinor 
  pbarAntiSpinor.SetP4(pbar4Vec, pbarMass);

  DiracSpinor pSpinor; // p spinor
  pSpinor.SetP4(p4Vec, pMass);

  DiracGamma5 gamma5;
  DiracGamma gamma; 

  std::vector< std::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  // fill singlet states
  
  Tensor<complex<double> > Psi0_m12_12=Bar(pbarAntiSpinor(1./2))*gamma5*pSpinor(-1./2);
  
  for(itJPCLS=jpclsSingletStates.begin(); itJPCLS!=jpclsSingletStates.end(); ++itJPCLS){
    OrbitalTensor currentTensor((*itJPCLS)->L);
    currentTensor.SetP4(pbar4Vec,p4Vec);

    PolVector currentPolVec((*itJPCLS)->J);
    currentPolVec.SetP4(allVec, allVec.M());

    Tensor<complex<double> > result= currentPolVec(0).Conjugate() | currentTensor;
    result*= Psi0_m12_12;

    _pbarpAmpSingletMap[(*itJPCLS)]=result(0);     
  }


  // fill triplet 0 states 
  Tensor<complex<double> > Psi1_m12_12=Bar(pbarAntiSpinor(-1./2))*gamma*pSpinor(1./2)-Bar(pbarAntiSpinor(-1./2))*cancelTerm*pSpinor(1./2);

  for(itJPCLS=jpclsTriplet0States.begin(); itJPCLS!=jpclsTriplet0States.end(); ++itJPCLS){
    OrbitalTensor currentTensor((*itJPCLS)->L);
    currentTensor.SetP4(pbar4Vec,p4Vec);

    PolVector currentPolVec((*itJPCLS)->J);
    currentPolVec.SetP4(allVec, allVec.M());

    Tensor<complex<double> > result;
    if ((*itJPCLS)->J > (*itJPCLS)->L) result= currentPolVec(0).Conjugate()  | (Psi1_m12_12 % currentTensor);
    else result= currentPolVec(0).Conjugate()  | (Psi1_m12_12 *  currentTensor );

    _pbarpAmpTriplet0Map[(*itJPCLS)]=result(0);  
  }

  // fill triplet 1 states
  LeviCivitaTensor lctTensor;
  Tensor<complex<double> > levipTensor=lctTensor*allVec;
  Tensor<complex<double> > Psi1_12_12=Bar(pbarAntiSpinor(1./2))*gamma*pSpinor(1./2)-Bar(pbarAntiSpinor(1./2))*cancelTerm*pSpinor(1./2);

  for(itJPCLS=jpclsTripletp1States.begin(); itJPCLS!=jpclsTripletp1States.end(); ++itJPCLS){
    OrbitalTensor currentTensor((*itJPCLS)->L);
    currentTensor.SetP4(pbar4Vec,p4Vec);

    PolVector currentPolVec((*itJPCLS)->J);
    currentPolVec.SetP4(allVec, allVec.M());

    Tensor<complex<double> > resultTr1;
    if ((*itJPCLS)->J > (*itJPCLS)->L){
      resultTr1= currentPolVec(1).Conjugate()  | (Psi1_12_12 % currentTensor);
    }
    else if ((*itJPCLS)->J == (*itJPCLS)->L){
      Tensor<complex<double> > tmpResultTr1=levipTensor*Psi1_12_12;
      Tensor<complex<double> > tmpResultTr1_1=tmpResultTr1.Permute(0, tmpResultTr1.Rank()-1)*currentTensor;
      resultTr1= currentPolVec(1).Conjugate()  | tmpResultTr1_1;

    }
    else{
      resultTr1= currentPolVec(1).Conjugate()  | (Psi1_12_12 *  currentTensor );
    }

    _pbarpAmpTripletp1Map[(*itJPCLS)]=resultTr1(0);

  }

  // fill triplet -1 states
  Tensor<complex<double> > Psi1_m12_m12=Bar(pbarAntiSpinor(-1./2))*gamma*pSpinor(-1./2)-Bar(pbarAntiSpinor(-1./2))*cancelTerm*pSpinor(-1./2);

  for(itJPCLS=jpclsTripletp1States.begin(); itJPCLS!=jpclsTripletp1States.end(); ++itJPCLS){
    OrbitalTensor currentTensor((*itJPCLS)->L);
    currentTensor.SetP4(pbar4Vec,p4Vec);

    PolVector currentPolVec((*itJPCLS)->J);
    currentPolVec.SetP4(allVec, allVec.M());

    Tensor<complex<double> > result;
    if ((*itJPCLS)->J > (*itJPCLS)->L){
      result= currentPolVec(-1).Conjugate()  | (Psi1_m12_m12 % currentTensor);
    }
    else if ((*itJPCLS)->J == (*itJPCLS)->L){
      Tensor<complex<double> > tmpResult=levipTensor*Psi1_m12_m12;
      Tensor<complex<double> > tmpResult_1=tmpResult.Permute(0, tmpResult.Rank()-1)*currentTensor; 
      result= currentPolVec(-1).Conjugate()  | tmpResult_1;      
    }
    else{
      result= currentPolVec(-1).Conjugate()  | (Psi1_m12_m12 *  currentTensor );
    }

    _pbarpAmpTripletm1Map[(*itJPCLS)]=result(0);
  } 
}


