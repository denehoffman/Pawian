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

#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/AbsEnv.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"

IsobarTensorDecay::IsobarTensorDecay(Particle* mother, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv) :
  IsobarLSDecay(mother, daughter1, daughter2, theEnv)
  ,_polMother(PolVector(_motherJPCPtr->J))
  ,_polDaughter1(PolVector(_daughter1JPCPtr->J))
  ,_polDaughter2(PolVector(_daughter2JPCPtr->J))
  ,_lctTensor(LeviCivitaTensor())
  , _metricTensor(MetricTensor())
{
}

IsobarTensorDecay::IsobarTensorDecay(std::shared_ptr<const jpcRes> motherJPCPtr, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv, std::string motherName) :
  IsobarLSDecay(motherJPCPtr, daughter1, daughter2, theEnv, motherName)
  ,_polMother(PolVector(_motherJPCPtr->J))
  ,_polDaughter1(PolVector(_daughter1JPCPtr->J))
  ,_polDaughter2(PolVector(_daughter2JPCPtr->J))
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

void IsobarTensorDecay::fillWignerDs(std::map<std::string, Vector4<double> >& fsMap, EvtData* evtData){
  int evtNo=evtData->evtNo;
  std::map<int, bool>::const_iterator it = _alreadyFilledMap.find(evtNo);
  if(it!=_alreadyFilledMap.end() &&  it->second) return; //already filled 

  if (!_daughter1IsStable) _absDecDaughter1->fillWignerDs(fsMap, evtData);
  if (!_daughter2IsStable) _absDecDaughter2->fillWignerDs(fsMap, evtData);

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
  //  DebugMsg << "daughter2_4Vec:\t"  << _daughter2->name() << "\t" << daughter2_4Vec << endmsg;
  
  Vector4<double> mother_4Vec=daughter1_4Vec+daughter2_4Vec;

  Spin spinMother=_motherJPCPtr->J;
  Spin spinDaughter1=_daughter1JPCPtr->J;
  Spin spinDaughter2=_daughter2JPCPtr->J;
  
  _polMother.SetP4(mother_4Vec,mother_4Vec.M());
  //  _polMother.SetP4(mother_4Vec,mother_4Vec.M());
  //  DebugMsg << "_polMother:\t" << _polMother << endmsg;   
  _polDaughter1.SetP4(daughter1_4Vec, daughter1_4Vec.M());
  _polDaughter2.SetP4(daughter2_4Vec, daughter2_4Vec.M());
  
  Spin lam12Max=spinDaughter1+spinDaughter2;
  if(lam12Max>spinMother) lam12Max=spinMother;
  
  Spin lamMotherMax=spinMother;
  if (!_hasMotherPart && spinMother>1) lamMotherMax=1; //attention: this is only valid for pbar p or e+ e- reactions; must be moved to individual specific classes

  std::vector< std::shared_ptr<const JPCLS> > theJPCLSAmps=JPCLSAmps();
  DebugMsg << name() << endmsg;  
  std::vector< std::shared_ptr<const JPCLS> >::iterator itJPCLS;
  for(itJPCLS=theJPCLSAmps.begin(); itJPCLS!=theJPCLSAmps.end(); ++itJPCLS){
    //      (*itJPCLS)->print(std::cout);
    //       std::cout << std::endl;
    Spin L=(*itJPCLS)->L;
    Spin S=(*itJPCLS)->S;
    int s1s2S=spinDaughter1+spinDaughter2+(*itJPCLS)->S;

    Tensor<complex<double> > leviPssTensor;
    bool add_lctForChi=true;
    if( s1s2S%2 ==0 ) { //odd
      add_lctForChi=false;
      leviPssTensor(0)=complex<double>(1.,0.);  
    }
    else{
      leviPssTensor=_lctTensor*mother_4Vec;
    }
    DebugMsg << "leviPssTensor: " << leviPssTensor << endmsg;

    int SLms1=(*itJPCLS)->S+(*itJPCLS)->L+_motherJPCPtr->J;

    Tensor<complex<double> > leviPlsTensor;
    bool add_lctForTensor=true;
    if( SLms1%2 ==0 ){ //odd
      add_lctForTensor=false;
      leviPlsTensor=complex<double>(1.,0.);  
    }
    else{
      leviPlsTensor=_lctTensor*mother_4Vec;
    }
    DebugMsg << "leviPlsTensor: " << leviPlsTensor << endmsg;

    OrbitalTensor orbTensor(L);
    orbTensor.SetP4(daughter1_4Vec, daughter2_4Vec);
    DebugMsg << "orbTensor:\t" << orbTensor << endmsg;    
    
    for (Spin lamMother=-lamMotherMax; lamMother<=lamMotherMax; ++lamMother){
      DebugMsg << "lamMother:\t" << lamMother << endmsg;
 
      Tensor<complex<double> > epsilonMotherProject = _polMother(lamMother);
      DebugMsg << "epsilonMotherProject:\t" << epsilonMotherProject << endmsg;

      // //calculate ls part;      
      // Tensor<complex<double> > lsPartTensor;
      // if(add_lctForTensor){
      // 	if(epsilonMotherProject.Rank() >= orbTensor.Rank()){
      // 	  lsPartTensor= epsilonMotherProject*leviPlsTensor;
      // 	  lsPartTensor.Permute(1,lsPartTensor.Rank());
      // 	  lsPartTensor=lsPartTensor|orbTensor;
      // 	}
      // 	else{
      // 	  lsPartTensor= orbTensor*leviPlsTensor;
      // 	  lsPartTensor.Permute(1,lsPartTensor.Rank());
      // 	  lsPartTensor=lsPartTensor|epsilonMotherProject;
      // 	}
      // }
      // for !add_lctForTensor calculation will be done in heli loops
 
      //calculate chi
      Tensor<complex<double> > chi12;
      Tensor<complex<double> > s12SpinProjector;
      PolVector part12PolVec(S);
      part12PolVec.SetP4(mother_4Vec,mother_4Vec.M());
      s12SpinProjector=part12PolVec.Projector();
      DebugMsg << "s12SpinProjector:\t" << s12SpinProjector << endmsg;


      for (Spin lamDaughter1=-spinDaughter1; lamDaughter1 <= spinDaughter1; ++lamDaughter1){
	DebugMsg << "lamDaughter1:\t" << lamDaughter1 << endmsg;

	Tensor<complex<double> > epsilonDaughter1Project = _polDaughter1(lamDaughter1);
	DebugMsg << "epsilonDaughter1Project:\t" << epsilonDaughter1Project << endmsg;

	for (Spin lamDaughter2=-spinDaughter2; lamDaughter2<=spinDaughter2; ++lamDaughter2){
	  DebugMsg << "lamDaughter2:\t" << lamDaughter2 << endmsg;
	  Tensor<complex<double> > epsilonDaughter2Project=_polDaughter2(lamDaughter2);
	  DebugMsg << "epsilonDaughter2Project:\t" << epsilonDaughter2Project << endmsg;

	  Tensor<complex<double> > chiPart;
	  if(add_lctForChi){
	    chiPart=(conj(epsilonDaughter1Project)*leviPssTensor)*conj(epsilonDaughter2Project);
	  }
	  else{
	    chiPart=conj(epsilonDaughter1Project)%conj(epsilonDaughter2Project);
	  } 
	  DebugMsg << "chiPart:\t" << chiPart << endmsg;

	  chi12=s12SpinProjector|chiPart;	  
	  DebugMsg << "chi12:\t" << chi12 << endmsg;

	  //calculate ls part;	  
	  Tensor<complex<double> > lsPartTensor;
	  if(add_lctForTensor){
	    if(epsilonMotherProject.Rank() >= orbTensor.Rank()){
	      lsPartTensor= epsilonMotherProject*leviPlsTensor;
	      lsPartTensor.Permute(1,lsPartTensor.Rank());
	      lsPartTensor=lsPartTensor|orbTensor;
	    }
	    else{
	      lsPartTensor= orbTensor*leviPlsTensor;
	      lsPartTensor.Permute(1,lsPartTensor.Rank());
	      lsPartTensor=lsPartTensor|epsilonMotherProject;
	    }
	  }
	  if(!add_lctForTensor){
	    if (chi12.Rank()==fabs(epsilonMotherProject.Rank()-orbTensor.Rank())) lsPartTensor= epsilonMotherProject|orbTensor; 
	    else if (chi12.Rank()==(epsilonMotherProject.Rank()+orbTensor.Rank())) lsPartTensor= epsilonMotherProject%orbTensor;
	    else if (epsilonMotherProject.Rank()==orbTensor.Rank()){
	      int noOfContractions=int((epsilonMotherProject.Rank()+orbTensor.Rank()-chi12.Rank())/2);
	      if (noOfContractions<=epsilonMotherProject.Rank() && noOfContractions>0) lsPartTensor= epsilonMotherProject.Contract(orbTensor, noOfContractions);
	    }
	    else{
	      Alert << "wrong ranks chi12.Rank()=" << chi12.Rank() 
		    << "\tepsilonMotherProject.Rank()=" << epsilonMotherProject.Rank()
		    << "\torbTensor.Rank()=" << orbTensor.Rank()
		    << endmsg;
	      exit(0);
	    }
	  }
	  DebugMsg << "lsPartTensor:\t" << lsPartTensor << endmsg;


	  Tensor<complex<double> > result=lsPartTensor|chi12;
	  DebugMsg << "result:\t" << result << endmsg;

	  if(result.Rank()>0){
	    Alert << "result.Rank()=" << result.Rank() << " > 0" << endmsg;
	    exit(0);
	  }
	  evtData->ComplexDouble5SpinString[_name][L][S][lamMother][lamDaughter1][lamDaughter2]=result(0);

	}
      }
    }
  }
  
}
