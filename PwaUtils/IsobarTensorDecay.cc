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

// IsobarTensorDecay class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>

#include "PwaUtils/AbsDynamics.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++Extension/PawianUtils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/DynRegistry.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "ConfigParser/ParserBase.hh"

IsobarTensorDecay::IsobarTensorDecay(Particle* mother, Particle* daughter1, Particle* daughter2, ChannelID channelID, std::string typeName) :
  IsobarLSDecay(mother, daughter1, daughter2, channelID, typeName)
  ,_polMother(PolVector(_motherIGJPCPtr->J))
  ,_polDaughter1(PolVector(_daughter1IGJPCPtr->J))
  ,_polDaughter2(PolVector(_daughter2IGJPCPtr->J))
  ,_lctTensor(LeviCivitaTensor())
  , _metricTensor(MetricTensor())
{
}

IsobarTensorDecay::IsobarTensorDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, Particle* daughter2, ChannelID channelID, std::string motherName, std::string typeName) :
  IsobarLSDecay(motherIGJPCPtr, daughter1, daughter2, channelID, motherName, typeName)
  ,_polMother(PolVector(_motherIGJPCPtr->J))
  ,_polDaughter1(PolVector(_daughter1IGJPCPtr->J))
  ,_polDaughter2(PolVector(_daughter2IGJPCPtr->J))
  ,_lctTensor(LeviCivitaTensor())
  ,_metricTensor(MetricTensor())
{
}

IsobarTensorDecay::~IsobarTensorDecay(){
}

void IsobarTensorDecay::print(std::ostream& os) const{
  os << "\nJPCLS tensor amplitudes for decay\t" << _name << ":\n";
  os << "suffix for fit parameter name:\t" << _fitParamSuffix << "\n";

  std::vector< std::shared_ptr<const JPCLS> >::const_iterator it;
  for (it = _JPCLSDecAmps.begin(); it!= _JPCLSDecAmps.end(); ++it){
    (*it)->print(os);
    os << "\n";
  }

  AbsDecay::print(os);
  os << "\n";
}

void IsobarTensorDecay::fillWignerDs(std::map<std::string, Vector4<double> >& fsMap, Vector4<double>& prodParticle4Vec, EvtData* evtData){
  int evtNo=evtData->evtNo;
  std::map<int, bool>::const_iterator it = _alreadyFilledMap.find(evtNo);
  if(it!=_alreadyFilledMap.end() &&  it->second) return; //already filled

  // if (!_daughter1IsStable) _absDecDaughter1->fillWignerDs(fsMap, evtData);
  // if (!_daughter2IsStable) _absDecDaughter2->fillWignerDs(fsMap, evtData);

  Vector4<double> all4Vec(0.,0.,0.,0.);
  //fill all4Vec
  std::map<std::string , Vector4<double> >::iterator itMap;
  for(itMap=fsMap.begin(); itMap!=fsMap.end(); ++itMap){
    all4Vec+=itMap->second;
  }

  //fill daughter1 and daughter2 4Vec
  Vector4<double> daughter1_4Vec(0.,0.,0.,0.);
  std::vector<Particle*> finalStatePart1;
  if(0!=_absDecDaughter1) finalStatePart1=_absDecDaughter1->finalStateParticles();
  else finalStatePart1.push_back(_daughter1);

  std::vector<Particle*>::iterator itParticle;

  for (itParticle= finalStatePart1.begin(); itParticle != finalStatePart1.end(); ++itParticle){
    itMap=fsMap.find( (*itParticle)->name() );
    daughter1_4Vec+=itMap->second;
  }

  //  DebugMsg << "daughter1_4Vec:\t" << _daughter1->name() << "\t" << daughter1_4Vec << endmsg;

  Vector4<double> daughter2_4Vec(0.,0.,0.,0.);

  std::vector<Particle*> finalStatePart2;
  if(0!=_absDecDaughter2) finalStatePart2=_absDecDaughter2->finalStateParticles();
  else finalStatePart2.push_back(_daughter2);

  for (itParticle= finalStatePart2.begin(); itParticle != finalStatePart2.end(); ++itParticle){
       itMap=fsMap.find( (*itParticle)->name() );
       daughter2_4Vec+=itMap->second;
  }

  if(_isProdAmp){
    if (!_daughter1IsStable) _absDecDaughter1->fillWignerDs(fsMap, daughter1_4Vec, evtData);
    if (!_daughter2IsStable) _absDecDaughter2->fillWignerDs(fsMap, daughter2_4Vec, evtData);
  }
  else{
    if (!_daughter1IsStable) _absDecDaughter1->fillWignerDs(fsMap, prodParticle4Vec, evtData);
    if (!_daughter2IsStable) _absDecDaughter2->fillWignerDs(fsMap, prodParticle4Vec, evtData);
  }


  //  DebugMsg << "daughter2_4Vec:\t"  << _daughter2->name() << "\t" << daughter2_4Vec << endmsg;

  Vector4<double> mother_4Vec=daughter1_4Vec+daughter2_4Vec;

  Vector4<double> daughter1Tensor4Vec=daughter1_4Vec;
  Vector4<double> daughter2Tensor4Vec=daughter2_4Vec;
  Vector4<double> motherTensor4Vec=mother_4Vec;

  //transorm into the CM frame of the initial two-body-system (e.g. pbar p, e+ e-, gamma p)
  daughter1Tensor4Vec.Boost(all4Vec);
  daughter2Tensor4Vec.Boost(all4Vec);
  motherTensor4Vec.Boost(all4Vec);

  Vector4<double> daughter2HelMother(0.,0.,0.,0.);
  // if(_hasMotherPart){
  //   if(mother_4Vec==prodParticle4Vec){
  //     daughter2HelMother=daughter2_4Vec;
  //     daughter2HelMother.Boost(daughter2HelMother); //is this correct????
  //   }
  //   else{
  //     daughter2HelMother=helicityVec(prodParticle4Vec, mother_4Vec, daughter2_4Vec);
  //   }
  // }
  // else{
  daughter2HelMother=daughter2_4Vec;    //only needed for prod amps
    daughter2HelMother.Boost(mother_4Vec);
    //  }


  //rotate everything into the cms flight direction 
  // Vector4<double> allRot4Vec=all4Vec;

  // allRot4Vec.RotateZ(-all4Vec.Phi());
  // allRot4Vec.RotateY(-all4Vec.Theta());
  // // // // std::cout << "\n\nall4Vec: " << all4Vec << std::endl;
  // // // // std::cout << "\nallRot4Vec: " << allRot4Vec << std::endl;

  // daughter1Tensor4Vec.RotateZ(-all4Vec.Phi());
  // daughter1Tensor4Vec.RotateY(-all4Vec.Theta());

  // daughter2Tensor4Vec.RotateZ(-all4Vec.Phi());
  // daughter2Tensor4Vec.RotateY(-all4Vec.Theta());

  // motherTensor4Vec.RotateZ(-all4Vec.Phi());
  // motherTensor4Vec.RotateY(-all4Vec.Theta());

  Spin spinMother=_motherIGJPCPtr->J;
  Spin spinDaughter1=_daughter1IGJPCPtr->J;
  Spin spinDaughter2=_daughter2IGJPCPtr->J;

  double motherMass=motherTensor4Vec.M();
  if (motherMass<1.e-5) motherMass=0.;
  else if(motherMass != motherMass) motherMass=0.;
  
  _polMother.SetP4(motherTensor4Vec, motherMass);


  double daughter1Mass=daughter1Tensor4Vec.M();
  double daughter2Mass=daughter2Tensor4Vec.M();
  if(daughter1Mass!=daughter1Mass || daughter1Mass < 1.e-5) daughter1Mass=0.;
  if(daughter2Mass!=daughter2Mass || daughter2Mass < 1.e-5) daughter2Mass=0.;

  _polDaughter1.SetP4(daughter1Tensor4Vec, daughter1Mass);
  _polDaughter2.SetP4(daughter2Tensor4Vec, daughter2Mass);

  Spin lamMotherMax=spinMother;
  if (!_hasMotherPart && spinMother>1){
    if ( GlobalEnv::instance()->Channel(_channelId)->channelType()==AbsChannelEnv::CHANNEL_PBARP 
         ||  GlobalEnv::instance()->Channel(_channelId)->channelType()==AbsChannelEnv::CHANNEL_EPEM ){
      lamMotherMax=1; //attention
    }
    else if (GlobalEnv::instance()->Channel(_channelId)->channelType()==AbsChannelEnv::CHANNEL_GAMMAP){
      lamMotherMax=3./2;
    }
  }

  std::vector< std::shared_ptr<const JPCLS> > theJPCLSAmps=JPCLSAmps();
  DebugMsg << name() << endmsg;
  std::vector< std::shared_ptr<const JPCLS> >::iterator itJPCLS;
  for(itJPCLS=theJPCLSAmps.begin(); itJPCLS!=theJPCLSAmps.end(); ++itJPCLS){
    //    (*itJPCLS)->print(std::cout);
    //   std::cout << std::endl;
    Spin L=(*itJPCLS)->L;
    Spin S=(*itJPCLS)->S;
    DebugMsg << "L: " << L << "\tS: " << S << endmsg;

    OrbitalTensor orbTensor(L);
    orbTensor.SetP4(daughter1Tensor4Vec, daughter2Tensor4Vec);
    DebugMsg << "orbTensor:\t" << orbTensor << endmsg;

    for (Spin lamMother=-lamMotherMax; lamMother<=lamMotherMax; ++lamMother){
      DebugMsg << "lamMother:\t" << lamMother << endmsg;

      Tensor<complex<double> > epsilonMotherProject = _polMother(lamMother);
      DebugMsg << "epsilonMotherProject:\t" << epsilonMotherProject << endmsg;

       //calculate chi
      Tensor<complex<double> > s12SpinProjector;
      PolVector part12PolVec(S);

      part12PolVec.SetP4(motherTensor4Vec, motherMass);
      s12SpinProjector=part12PolVec.Projector();
      DebugMsg << "s12SpinProjector:\t" << "should have rank " << 2*(*itJPCLS)->S <<"\trank: " << s12SpinProjector.Rank() << endmsg;
      if(s12SpinProjector.Rank()<=2) DebugMsg << s12SpinProjector << endmsg;


      for (Spin lamDaughter1=-spinDaughter1; lamDaughter1 <= spinDaughter1; ++lamDaughter1){
	DebugMsg << "lamDaughter1:\t" << lamDaughter1 << endmsg;

	Tensor<complex<double> > epsilonDaughter1Project = _polDaughter1(lamDaughter1);
	if(_isDaughter1Photon){
	  epsilonDaughter1Project = epsilonDaughter1Project
	    - (daughter2Tensor4Vec*epsilonDaughter1Project)/(daughter1Tensor4Vec*daughter2Tensor4Vec)*daughter1Tensor4Vec;
	}
	DebugMsg << "epsilonDaughter1Project:\t" << epsilonDaughter1Project << endmsg;

	for (Spin lamDaughter2=-spinDaughter2; lamDaughter2<=spinDaughter2; ++lamDaughter2){
	  // Spin lambda = lamDaughter1-lamDaughter2;
	  // if( fabs(lambda)>(*itJPCLS)->J || fabs(lambda)>(*itJPCLS)->S) continue;
 	  
	  DebugMsg << "lamDaughter2:\t" << lamDaughter2 << endmsg; 

	  Tensor<complex<double> > epsilonDaughter2Project=_polDaughter2(lamDaughter2);
	  if(_isDaughter2Photon){
	    epsilonDaughter2Project = epsilonDaughter2Project
	      - (daughter1Tensor4Vec*epsilonDaughter2Project)/(daughter2Tensor4Vec*daughter1Tensor4Vec)*daughter2Tensor4Vec;
	  }

	  DebugMsg << "epsilonDaughter2Project:\t" << epsilonDaughter2Project << endmsg;


	  Tensor<complex<double> > chi12;
	  calcChi12(s12SpinProjector, epsilonDaughter1Project, epsilonDaughter2Project, motherTensor4Vec, chi12);

	  //calculate ls part;
	  Tensor<complex<double> > lsPartTensorNew;
	  calcLSpart(orbTensor, chi12, S, motherTensor4Vec, lsPartTensorNew);

	  DebugMsg << "epsilonDaughter2Project rank:\t" << epsilonDaughter2Project.Rank() << endmsg;
	  if(epsilonDaughter2Project.Rank()<3) DebugMsg << epsilonDaughter2Project << endmsg;

	  DebugMsg << "lsPartTensorNew rank:\t" << lsPartTensorNew.Rank() << endmsg;
	  if(lsPartTensorNew.Rank()<3) DebugMsg << lsPartTensorNew << endmsg;
 
	  Tensor<complex<double> > result=epsilonMotherProject|lsPartTensorNew;
	  DebugMsg << "result " << name()<<": " << result << "\n" << endmsg;

	  if(result.Rank()>0){
	    Alert << "result.Rank()=" << result.Rank() << " > 0" << endmsg;
	    Alert << "decay name: " <<  name() << endmsg;
	    exit(0);
	  }
	  if (norm(result(0)) >1000.){
	    DebugMsg << "norm(result(0)) (L="<< L << ", S=" << S << ", lamMother=" << lamMother << ", lamDaughter1=" << lamDaughter1
		     << ", lamDaughter2=" << lamDaughter2  << ")= " << norm(result(0)) << endmsg;
	  }
	  Id3StringType IdLamXLam1Lam2=FunctionUtils::spin3Index(lamMother, lamDaughter1, lamDaughter2);
	  unsigned short currentLSId=(*itJPCLS)->idnumberLS;
	  evtData->ComplexLS3Spin[_nameId][currentLSId][IdLamXLam1Lam2]=result(0);	  
	}
      }
    }
  }

  bool fillqVals=false;
  if(_isProdAmp && _useProdBarrier) fillqVals=true;
  else if(0!=_absDynPtr){
    if(_absDynPtr->type()=="BlattWBarrierTensorDynamics") fillqVals=true; 
  }
  
  if(fillqVals){
    double qVal=daughter2HelMother.P();
    double qValNorm=PawianQFT::breakupMomQDefault(mother_4Vec.M(), massSumFsParticlesDec1(), massSumFsParticlesDec2()).real();
    evtData->DoubleMassId[_wignerDqId]=qVal;
    evtData->DoubleMassId[_wignerDqNormId] = qValNorm;
  }
  _alreadyFilledMap[evtNo]=true;  
}

void IsobarTensorDecay::calcChi12(Tensor<complex<double> >& s12SpinProjector, Tensor<complex<double> >& epsilonDaughter1Project, Tensor<complex<double> >& epsilonDaughter2Project, Vector4<double>& mother_4Vec, Tensor<complex<double> >& result){
  Spin spin12=s12SpinProjector.Rank()/2.;
  Spin spin1=epsilonDaughter1Project.Rank();
  Spin spin2=epsilonDaughter2Project.Rank();

  int s1s2S=spin1+spin2+spin12;

  Tensor<complex<double> > leviPssTensor;
  bool add_lctForChi=true;
  if( s1s2S%2 ==0 ) { //even
    add_lctForChi=false;
    leviPssTensor(0)=complex<double>(1.,0.);
  }
  else{
    leviPssTensor=_lctTensor*mother_4Vec;
  }
  DebugMsg << "leviPssTensor: " << leviPssTensor << endmsg;

  Tensor<complex<double> > chiPart;
  if(add_lctForChi){
    chiPart=leviPssTensor*conj(epsilonDaughter1Project);
    chiPart.Permute(0,chiPart.Rank()-1);
    //    chiPart=conj(epsilonDaughter1Project)*leviPssTensor;
    chiPart=chiPart*conj(epsilonDaughter2Project);
  }
  else{
    if(spin12==0) chiPart=conj(epsilonDaughter1Project)|conj(epsilonDaughter2Project); 
    else chiPart=conj(epsilonDaughter1Project)%conj(epsilonDaughter2Project);
  }
  DebugMsg << "chiPart:\t" << chiPart << endmsg;

  //  result=s12SpinProjector|chiPart;

  int noContractionsChi12=0.5*(spin12+spin1+spin2);
  int noRemainingContr=noContractionsChi12;
  if (add_lctForChi){
    noContractionsChi12=0.5*(spin12+4+1+spin1+spin2);
    noRemainingContr=noContractionsChi12-3;    
  }

  result=s12SpinProjector.Contract(chiPart, noRemainingContr);

  DebugMsg << "no of contractions chi12:\t" << noContractionsChi12 << endmsg;
  int rankCh12=2*spin12+spin1+spin2-2*noContractionsChi12;
  if(add_lctForChi) rankCh12=2*spin12+spin1+spin2+5-2*noContractionsChi12; 
  DebugMsg << "chi12:\t" << "should have rank " << rankCh12 << " and (cross check): " << spin12 <<"\trank" << result.Rank() << endmsg;
  if(result.Rank()<3){
    DebugMsg << result << endmsg;
  }
}


void IsobarTensorDecay::calcLSpart(OrbitalTensor& orbTensor, Tensor<complex<double> >& chi12, Spin spin12, Vector4<double>& mother_4Vec, Tensor<complex<double> >& result){
  int orbitalL=orbTensor.Rank();
  if(orbTensor.Rank()<3){
    DebugMsg << "orbTensor: " << orbTensor << endmsg;    
  }


  int s12LJ=spin12+orbitalL+_motherIGJPCPtr->J; 
  Tensor<complex<double> > leviPlsTensor;
  bool add_lct=true;
  int noOfContractions=0;
  if( s12LJ%2 ==0 ){ //even
    add_lct=false;
    leviPlsTensor=complex<double>(1.,0.);
    noOfContractions=(orbitalL+spin12-_motherIGJPCPtr->J)/2;
  }
  else{
    leviPlsTensor=_lctTensor*mother_4Vec;
    noOfContractions=(4+1+orbitalL+spin12-_motherIGJPCPtr->J)/2;
  }
  DebugMsg << "leviPlsTensor: " << leviPlsTensor << endmsg;
  DebugMsg << "no of contr.: " << noOfContractions << endmsg;

  Tensor<complex<double> > lsPartTensor;
  if(add_lct){
    lsPartTensor=leviPlsTensor*orbTensor;
    lsPartTensor.Permute(0,lsPartTensor.Rank()-1);

    //    lsPartTensor=orbTensor*leviPlsTensor;

    //     DebugMsg << "lsPartTensor w levi and permute: rank" << lsPartTensor.Rank() << endmsg;
    //      DebugMsg << "chi12.Rank() " << chi12.Rank() << endmsg;
    result=lsPartTensor.Contract(chi12, noOfContractions-2);
  }
  else{ // add_lct=false
    if(noOfContractions==0) result=orbTensor%chi12;
    else if(orbTensor.Rank()>=chi12.Rank()){
      result=orbTensor.Contract(chi12, noOfContractions); 
    }
    else{
      result=chi12.Contract(orbTensor, noOfContractions);
    } 
  }

    DebugMsg << "LS part: should have rank " << _motherIGJPCPtr->J << "\t rank" << result.Rank() << endmsg;
    if(result.Rank()<3){
      DebugMsg << result << endmsg;
    }
}

void IsobarTensorDecay::enableProdBarrier(){
  if(_dynEnabled){
    Alert << "dynamics already enabled for " << name() << endmsg;
    exit(1);
  }
  if(!_prodChannelInfo->isProductionChannel()){
    WarningMsg << name() << " is not a production amplitide! Barrier factors for the production can not be enabled!" << endmsg;
    return;
  }
  if(!_prodChannelInfo->withProdBarrier()){
    WarningMsg << name() << "production barrier disabled" << endmsg;
    return;
  }
  if(_prodChannelInfo->prodBarrierType()!="BlattWBarrierTensor"){
    Alert << name() << "production barrier enabled with type " << _prodChannelInfo->prodBarrierType() 
	  << "\tIsobarTensorDecay only supports type BlattWBarrierTensor" << endmsg;
    exit(1);
  }
  _useProdBarrier=true;
  _dynType="BlattWBarrierTensor";
  _qR=_prodChannelInfo->qRPod();
  InfoMsg << "Barrier factors for production amplitude " << name() << " enabled!" << endmsg;
  _absDynPtr=DynRegistry::instance()->getDynamics(shared_from_this()); 
  _dynEnabled=true;
}
