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

// IsobarTensorPsiToGamXDecay class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>

#include "PwaUtils/IsobarTensorPsiToGamXDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsDynamics.hh"
#include "PwaUtils/DynRegistry.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "ConfigParser/ParserBase.hh"

IsobarTensorPsiToGamXDecay::IsobarTensorPsiToGamXDecay(Particle* mother, Particle* daughter1_gamma, Particle* daughter2, ChannelID channelID) :
  IsobarTensorDecay(mother, daughter1_gamma, daughter2, channelID, "IsobarTensorPsiToGamXDecay")
  ,_XisEven(false)
  ,_noOfAmps(0)
{
  if(mother->twoIso() !=0 || mother->theParity() != -1 || mother->twoJ() != 2 || mother->theParity() != -1 ||  mother->theCParity() || mother->mass() < 0.01){
    Alert << "mother particle does not have the quantum number combinations JPC=1-- or its mass is below 0.01 GeV/c2 !!!" << endmsg;
    exit(0); 
  } 

  if(daughter1_gamma->twoJ() != 2 || daughter1_gamma->theParity() != -1 ||  daughter1_gamma->theCParity() || daughter1_gamma->mass() > 1.e-6){
    Alert << "daughter particle 1 is not a photon" << endmsg;
    exit(0); 
  } 
  if(daughter2->twoJ()%4 == 0) _XisEven=true;

  if(!_XisEven && mother->theParity() == -1 && daughter2->twoJ()>=6){
    Alert << "radiative decays of Psi -> odd^- + gamma with J>=3 are not supported so far" << endmsg;
    exit(0);
  }

  if(!_XisEven && mother->theParity() == 1 && daughter2->twoJ()>=6){
    Alert << "radiative decays of Psi -> odd^+ + gamma with J>=3 are not supported so far" << endmsg;
    exit(0);
  }

  if (daughter2->twoJ()==0) _noOfAmps=1;
  else if(daughter2->twoJ()==2) _noOfAmps=2;
  else _noOfAmps=3;  

  Info << "daughter2->twoJ()%2 = " << daughter2->twoJ()%2 << "_XisEven: " << _XisEven <<  "  daughter2->theParity(): " << daughter2->theParity() << endmsg;
  fillAmpLMap();  
}

IsobarTensorPsiToGamXDecay::IsobarTensorPsiToGamXDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1_gamma, Particle* daughter2, ChannelID channelID, std::string motherName) :
  IsobarTensorDecay(motherIGJPCPtr, daughter1_gamma, daughter2, channelID, motherName, "IsobarTensorPsiToGamXDecay")
  ,_XisEven(false)
{ 
  if(motherIGJPCPtr->I != 0 || motherIGJPCPtr->G != -1 ||  motherIGJPCPtr->J != 1 || motherIGJPCPtr->P != -1 || motherIGJPCPtr->C != -1 ){
    Alert << "mother particle does not have the quantum number combinations IG(JPC)=0-(1--) !!!" << endmsg;
    exit(0);
  }

  if(daughter1_gamma->twoJ() != 2 || daughter1_gamma->theParity() != -1 ||  daughter1_gamma->theCParity() !=-1 || daughter1_gamma->mass() > 1.e-6){
    Alert << "daughter particle 1 is not a photon" << endmsg;
    Alert << "daughter1_gamma->twoJ(): " << daughter1_gamma->twoJ() << "\n"
	  << "daughter1_gamma->theParity(): " << daughter1_gamma->theParity() << "\n"
	  << "daughter1_gamma->theCParity(): " << daughter1_gamma->theCParity() << "\n"
	  << "daughter1_gamma->mass(): " << daughter1_gamma->mass() << endmsg;
    exit(0); 
  }

  if(daughter2->twoJ()%4 ==0) _XisEven=true;

  if(!_XisEven && motherIGJPCPtr->P == -1 && daughter2->twoJ()>=6){
    Alert << "radiative decays of Psi -> odd^- + gamma with J>=3 are not supported so far" << endmsg;
    exit(0);
  }

  if(!_XisEven && motherIGJPCPtr->P == 1 && daughter2->twoJ()>=6){
    Alert << "radiative decays of Psi -> odd^+ + gamma with J>=3 are not supported so far" << endmsg;
    exit(0);
  }

  if (daughter2->twoJ()==0) _noOfAmps=1;
  else if(daughter2->twoJ()==2) _noOfAmps=2;
  else _noOfAmps=3;

  Info << "daughter2->twoJ()%2 = " << daughter2->twoJ()%2 << "_XisEven: " << _XisEven <<  "  daughter2->theParity(): " << daughter2->theParity() << endmsg;
  fillAmpLMap();  
}

IsobarTensorPsiToGamXDecay::~IsobarTensorPsiToGamXDecay(){
}

void IsobarTensorPsiToGamXDecay::print(std::ostream& os) const{
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

void IsobarTensorPsiToGamXDecay::fillWignerDs(std::map<std::string, Vector4<double> >& fsMap, Vector4<double>& prodParticle4Vec, EvtData* evtData){

 int evtNo=evtData->evtNo;
  std::map<int, bool>::const_iterator it = _alreadyFilledMap.find(evtNo);
  if(it!=_alreadyFilledMap.end() &&  it->second) return; //already filled

  Vector4<double> all4Vec(0.,0.,0.,0.);
  //fill all4Vec
  std::map<std::string , Vector4<double> >::iterator itMap;
  for(itMap=fsMap.begin(); itMap!=fsMap.end(); ++itMap){
    all4Vec+=itMap->second;
  }

  //fill daughter1 and daughter2 4Vec
  Vector4<double> daughter1Gam_4Vec(0.,0.,0.,0.);

  itMap=fsMap.find( _daughter1->name() );
  if(itMap!=fsMap.end())  daughter1Gam_4Vec=itMap->second;
  else {
    Alert << "particle with the name " << _daughter1->name() << " can not be found in the list of final state particles!!!" << endmsg;
    exit(0);     
  }

  Vector4<double> daughter2_4Vec(0.,0.,0.,0.);

  std::vector<Particle*> finalStatePart2;
  if(0!=_absDecDaughter2) finalStatePart2=_absDecDaughter2->finalStateParticles();
  else finalStatePart2.push_back(_daughter2);

  std::vector<Particle*>::iterator itParticle;
  for (itParticle= finalStatePart2.begin(); itParticle != finalStatePart2.end(); ++itParticle){
       itMap=fsMap.find( (*itParticle)->name() );
       daughter2_4Vec+=itMap->second;
  }

  if(_isProdAmp){
    if (!_daughter2IsStable) _absDecDaughter2->fillWignerDs(fsMap, daughter2_4Vec, evtData);
  }
  else{
    if (!_daughter2IsStable) _absDecDaughter2->fillWignerDs(fsMap, prodParticle4Vec, evtData);
  }


  Vector4<double> mother_4Vec=daughter1Gam_4Vec+daughter2_4Vec;

  Vector4<double> daughter1GamTensor4Vec=daughter1Gam_4Vec;
  Vector4<double> daughter2Tensor4Vec=daughter2_4Vec;
  Vector4<double> motherTensor4Vec=mother_4Vec;

  daughter1GamTensor4Vec.Boost(all4Vec);
  daughter2Tensor4Vec.Boost(all4Vec);
  motherTensor4Vec.Boost(all4Vec);

  Vector4<double> daughter2HelMother(0.,0.,0.,0.);
  if(_hasMotherPart){
    if(fabs(mother_4Vec==all4Vec)){
      daughter2HelMother=daughter2_4Vec;
      daughter2HelMother.Boost(daughter2HelMother); //is this correct????
    }
    else daughter2HelMother=helicityVec(prodParticle4Vec, mother_4Vec, daughter2_4Vec);
  }
  else{
    daughter2HelMother=daughter2_4Vec;
    daughter2HelMother.Boost(mother_4Vec);
  }

  Spin spinMother(1);
  Spin spinDaughter1Gam(1);
  Spin spinDaughter2=_daughter2IGJPCPtr->J;

  double motherMass=motherTensor4Vec.M();
  if (motherMass<1.e-5) motherMass=0.;
  else if(motherMass != motherMass) motherMass=0.;
  
  _polMother.SetP4(motherTensor4Vec, motherMass);


  double daughter1GamMass=0.;
  double daughter2Mass=daughter2Tensor4Vec.M();
  if(daughter2Mass!=daughter2Mass || daughter2Mass < 1.e-5) daughter2Mass=0.;

  _polDaughter1.SetP4(daughter1GamTensor4Vec, daughter1GamMass);
  _polDaughter2.SetP4(daughter2Tensor4Vec, daughter2Mass);
  Tensor<complex<double> > S_alpha_beta= (_lctTensor*daughter1GamTensor4Vec)*motherTensor4Vec;

  Tensor<complex<double> > currentAmp1;
  Tensor<complex<double> > currentAmp2;
  Tensor<complex<double> > currentAmp3;

  Spin lamMotherMax=spinMother;

  for (Spin lamMother=-lamMotherMax; lamMother<=lamMotherMax; ++lamMother){
    for (Spin lam1Gam=-1; lam1Gam<=1; lam1Gam+=2){
      Tensor<complex<double> > PsiGamPolProj=_polMother(lamMother)%conj(_polDaughter1(lam1Gam));

      for (Spin mX=-spinDaughter2; mX<=spinDaughter2; ++mX){
	Tensor<complex<double> > PpsySpinProjX_jg0=conj(_polDaughter2(mX));
	Tensor<complex<double> > PpsySpinProjX_jg1=conj(_polDaughter2(mX));	
	Tensor<complex<double> > PpsySpinProjX_jg2=conj(_polDaughter2(mX));

	if(spinDaughter2>0){
	  //J contractions
	  for (int i=0; i<spinDaughter2; ++i) PpsySpinProjX_jg0*=motherTensor4Vec;
	}
	if(spinDaughter2>1){
	  //J-1 contractions
	  for (int i=0; i<spinDaughter2-1; ++i) PpsySpinProjX_jg1*=motherTensor4Vec;
	}

	if(spinDaughter2>2){
	  //J-2 contractions
	  for (int i=0; i<spinDaughter2-2; ++i) PpsySpinProjX_jg2*=motherTensor4Vec;
	}

	DebugMsg << "PpsySpinProjX_jg0.Rank(): " << PpsySpinProjX_jg0.Rank() << "\nval: " << PpsySpinProjX_jg0 << endmsg;
	DebugMsg << "PpsySpinProjX_jg1.Rank(): " << PpsySpinProjX_jg1.Rank() << "\nval: " << PpsySpinProjX_jg1 << endmsg;
	DebugMsg << "PpsySpinProjX_jg2.Rank(): " << PpsySpinProjX_jg2.Rank() << "\nval: " << PpsySpinProjX_jg2 << endmsg;

	//case J^PC(X)=0+, 2+, 4+, ... and 1-, 3-, ...
	if((_XisEven && _daughter2IGJPCPtr->P==1) || (!_XisEven && _daughter2IGJPCPtr->P==-1) ){ //0++, 2++, 4++, ...
	  //first amplitude
	  MetricTensor g_munu;
 	  
	  currentAmp1=(PsiGamPolProj | g_munu) * PpsySpinProjX_jg0;
	  //	  if(_daughter2->twoJ()>=4) currentAmp1 *= PpsySpinProjX_jg2 | (motherTensor4Vec%motherTensor4Vec);

	  DebugMsg << "currentAmp1.Rank(): " << currentAmp1.Rank() << "\nval: " << currentAmp1 << endmsg;

	  Id3StringType IdLamMotherLamGamLamX=FunctionUtils::spin3Index(lamMother, lam1Gam, mX);
	  //	  evtData->ComplexDoubleInt3SpinString[_name][0][lamMother][lam1Gam][mX]=currentAmp1(0);	  
	  evtData->ComplexN3Spin[_nameId][0][IdLamMotherLamGamLamX]=currentAmp1(0);
	  if(_noOfAmps>1){
	    //second amplitude
	    currentAmp2=(_polMother(lamMother)*daughter1GamTensor4Vec)*( conj(_polDaughter1(lam1Gam)) | PpsySpinProjX_jg1 ); 

	    DebugMsg << "currentAmp2.Rank(): " << currentAmp2.Rank() << "\nval: " << currentAmp2 << endmsg;
	    
	    //evtData->ComplexDoubleInt3SpinString[_name][1][lamMother][lam1Gam][mX]=currentAmp2(0);	    
	    evtData->ComplexN3Spin[_nameId][1][IdLamMotherLamGamLamX]=currentAmp2(0);
	    if(_noOfAmps>2){
	      //3. amplitude
	      currentAmp3= PsiGamPolProj | PpsySpinProjX_jg2;
	      DebugMsg << "currentAmp3.Rank(): " << currentAmp3.Rank() << "\nval: " << currentAmp3 << endmsg;
	      //	      evtData->ComplexDoubleInt3SpinString[_name][2][lamMother][lam1Gam][mX]=currentAmp3(0);
	      evtData->ComplexN3Spin[_nameId][2][IdLamMotherLamGamLamX]=currentAmp3(0);
	    }
	  }  	  
	}

	// //case J^PC(X)=1+, 3+, 5+, ...
	// if(!_XisEven && _daughter2IGJPCPtr->P==1){
	    
	//     //1. amplitude
	//     Tensor<complex<double> > U_numu1= (_lctTensor*PpsySpinProjX_jg1)*motherTensor4Vec;
	//     currentAmp1 = PsiGamPolProj | U_numu1;
	//     DebugMsg << "currentAmp1.Rank(): " << currentAmp1.Rank() << "\nval: " << currentAmp1(0) << endmsg;
	//     evtData->ComplexDoubleInt3SpinString[_name][0][lamMother][lam1Gam][mX]=currentAmp1(0);

	//     //2. amplitude
	//     Tensor<complex<double> > U_numu2= daughter1GamTensor4Vec % (S_alpha_beta*PpsySpinProjX_jg1);
	//     currentAmp2=PsiGamPolProj | U_numu2;
	//     DebugMsg << "currentAmp2.Rank(): " << currentAmp2.Rank() << "\nval: " << currentAmp2(0) << endmsg;
	//     evtData->ComplexDoubleInt3SpinString[_name][1][lamMother][lam1Gam][mX]=currentAmp2(0);

	//     if(spinDaughter2>1){
	//       //3. amplitude ?????
	//       currentAmp2=complex<double>(0.,0.);
	//       evtData->ComplexDoubleInt3SpinString[_name][2][lamMother][lam1Gam][mX]=complex<double>(0.,0.);
	//     }
	// }

	//case J^PC(X)=0-, 2-, 4-, ... and 1+, 3+, ...
	if((_XisEven && _daughter2IGJPCPtr->P==-1) || (!_XisEven && _daughter2IGJPCPtr->P==1)){ //0-, 2-, 4-, ...
	  //if(_XisEven && _daughter2IGJPCPtr->P==-1){
	    //1. amplitude 
	    Tensor<complex<double> > U_numu1;
            U_numu1 = S_alpha_beta*PpsySpinProjX_jg0;
	    currentAmp1=PsiGamPolProj | U_numu1;
	    DebugMsg << "currentAmp1.Rank(): " << currentAmp1.Rank() << "\nval: " << currentAmp1(0) << endmsg;
	    //	    evtData->ComplexDoubleInt3SpinString[_name][0][lamMother][lam1Gam][mX]=currentAmp1(0);

	    Id3StringType IdLamMotherLamGamLamX=FunctionUtils::spin3Index(lamMother, lam1Gam, mX);
	  //	  evtData->ComplexDoubleInt3SpinString[_name][0][lamMother][lam1Gam][mX]=currentAmp1(0);	  
	    evtData->ComplexN3Spin[_nameId][0][IdLamMotherLamGamLamX]=currentAmp1(0);
	    
	    if(spinDaughter2>0){
	      //2. amplitude 
	      Tensor<complex<double> > U_numu2= daughter1GamTensor4Vec % (S_alpha_beta*PpsySpinProjX_jg1);
	      currentAmp2=PsiGamPolProj | U_numu2;
	      DebugMsg << "currentAmp2.Rank(): " << currentAmp2.Rank() << "\nval: " << currentAmp2(0) << endmsg;
	      //	      evtData->ComplexDoubleInt3SpinString[_name][1][lamMother][lam1Gam][mX]=currentAmp2(0);
	      evtData->ComplexN3Spin[_nameId][1][IdLamMotherLamGamLamX]=currentAmp2(0);

	      if(spinDaughter2>1){	      
		//3. amplitude ?????
		Tensor<complex<double> > U_numu3= (_lctTensor*(PpsySpinProjX_jg2*daughter1GamTensor4Vec))*motherTensor4Vec;
		//		Tensor<complex<double> > U_numu3= _lctTensor | PpsySpinProjX_jg2;
		currentAmp3=PsiGamPolProj | U_numu3;
		//		evtData->ComplexDoubleInt3SpinString[_name][2][lamMother][lam1Gam][mX]=currentAmp3(0);
		evtData->ComplexN3Spin[_nameId][2][IdLamMotherLamGamLamX]=currentAmp3(0);
	      }
	    }
	}

	// //case J^PC(X)=1-, 3-, 5-, ...
	// if(!_XisEven && _daughter2IGJPCPtr->P==-1){ //??????????????
	//   MetricTensor g_munu;
	  
	//   //1. amplitude
	//   DebugMsg << "S_alpha_beta*PpsySpinProjX_jg1: " << S_alpha_beta*PpsySpinProjX_jg1 << endmsg;
	//   DebugMsg << "(S_alpha_beta*PpsySpinProjX_jg1)*motherTensor4Vec: " << (S_alpha_beta*PpsySpinProjX_jg1)*motherTensor4Vec << endmsg;
	//   DebugMsg << "S_alpha_beta*motherTensor4Vec: " << S_alpha_beta*motherTensor4Vec << endmsg;
	//   DebugMsg << "(S_alpha_beta*motherTensor4Vec)*PpsySpinProjX_jg1: " << (S_alpha_beta*motherTensor4Vec)*PpsySpinProjX_jg1 << endmsg;
	//   DebugMsg << "(S_alpha_beta*PpsySpinProjX_jg1)*daughter1GamTensor4Vec: " << (S_alpha_beta*PpsySpinProjX_jg1)*daughter1GamTensor4Vec << endmsg;
	//   Tensor<complex<double> > U_numu1= g_munu*((S_alpha_beta*PpsySpinProjX_jg1)*motherTensor4Vec);
	//   DebugMsg << "U_numu1: " << U_numu1 << endmsg;

	//   currentAmp1 = PsiGamPolProj | U_numu1;
	//   DebugMsg << "currentAmp1.Rank(): " << currentAmp1.Rank() << "\nval: " << currentAmp1(0) << endmsg;
	//   evtData->ComplexDoubleInt3SpinString[_name][0][lamMother][lam1Gam][mX]=currentAmp1(0);

	//   //2. amplitude
	//   DebugMsg << "((_lctTensor*PpsySpinProjX_jg1)*motherTensor4Vec)*motherTensor4Vec: " << ((_lctTensor*PpsySpinProjX_jg1)*motherTensor4Vec)*motherTensor4Vec << endmsg;
	//   DebugMsg << "(_lctTensor*PpsySpinProjX_jg1)*motherTensor: " << (_lctTensor*PpsySpinProjX_jg1)*motherTensor4Vec << endmsg;

	//   Tensor<complex<double> > U_numu2= daughter1GamTensor4Vec % (((_lctTensor*PpsySpinProjX_jg1)*motherTensor4Vec)*motherTensor4Vec);
	//   currentAmp2=PsiGamPolProj | U_numu2;
	//   DebugMsg << "currentAmp2.Rank(): " << currentAmp2.Rank() << "\nval: " << currentAmp2(0) << endmsg;
	//   evtData->ComplexDoubleInt3SpinString[_name][1][lamMother][lam1Gam][mX]=currentAmp2(0);
	  
	//   if(spinDaughter2>1){
	//     //3. amplitude ?????
	//     currentAmp2=complex<double>(0.,0.);
	//     evtData->ComplexDoubleInt3SpinString[_name][2][lamMother][lam1Gam][mX]=complex<double>(0.,0.);
	//   }
	// }
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
    double qValNorm=breakupMomQ(mother_4Vec.M(), massSumFsParticlesDec1(), massSumFsParticlesDec2()).real();
    evtData->DoubleMassId[_wignerDqId]=qVal;
    evtData->DoubleMassId[_wignerDqNormId]=qValNorm;
  }
}


void IsobarTensorPsiToGamXDecay::fillAmpLMap(){

//case J^PC(X)=0+, 2+, 4+, ... and 1-, 3-, ...
  if((_XisEven && _daughter2IGJPCPtr->P==1) || (!_XisEven && _daughter2IGJPCPtr->P==-1) ){
    _ampLMap[0]=_daughter2IGJPCPtr->J;
    if(_noOfAmps>1) _ampLMap[1]=_daughter2IGJPCPtr->J;
    if(_noOfAmps>2) _ampLMap[2]=_daughter2IGJPCPtr->J-2; 
  }
  else{	//case J^PC(X)=0-, 2-, 4-, ... and 1+, 3+, ...
    _ampLMap[0]=_daughter2IGJPCPtr->J+1;
    if(_noOfAmps>1) _ampLMap[1]=_daughter2IGJPCPtr->J+1;
    if(_noOfAmps>2) _ampLMap[2]=_daughter2IGJPCPtr->J-1;
  }
}

void IsobarTensorPsiToGamXDecay::enableProdBarrier(){
  if(_dynEnabled){
    Alert << "dynamics already enabled for " << name() << endmsg;
    exit(1);
  }
  if(!_prodChannelInfo->isProductionChannel()){
    Warning << name() << " is not a production amplitide! Barrier factors for the production can not be enabled!" << endmsg;
    return;
  }
  if(!_prodChannelInfo->withProdBarrier()){
    Warning << name() << "production barrier disabled" << endmsg;
    return;
  }
  if(_prodChannelInfo->prodBarrierType()!="BlattWBarrierTensor"){
    Alert << name() << "production barrier enabled with type " << _prodChannelInfo->prodBarrierType() 
	  << "\tIsobarTensorPsiToGamXDecay only supports type BlattWBarrierTensor" << endmsg;
    exit(1);
  }


  if(!_isProdAmp){
    Alert << name() << " is not a production amplitide! Barrier factors for the production can not be enabled!" << endmsg;
    exit(1);
  }
  _useProdBarrier=true;
  _dynType="BlattWBarrierTensor";
  _qR=_prodChannelInfo->qRPod();
  Info << "Barrier factors for production amplitude " << name() << " enabled!" << endmsg;
  _absDynPtr=DynRegistry::instance()->getDynamics(shared_from_this()); 
  _dynEnabled=true;
}

