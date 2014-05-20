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

#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"

#include "Examples/Tutorial/DfuncClebschG/DecAngularDist.hh"
#include "PwaUtils/DataUtils.hh"
#include "qft++Extension/SphericalHarmonic.hh"


#include "ErrLogger/ErrLogger.hh"

DecAngularDist::DecAngularDist(std::shared_ptr<const jpcRes> motherJPC, std::shared_ptr<const jpcRes> daughter1JPC, std::shared_ptr<const jpcRes> daughter2JPC, bool sameDaughters) :
  _motherJPC(motherJPC)
  ,_daughter1JPC(daughter1JPC)
  ,_daughter2JPC(daughter2JPC)
  ,_sameDaughters(sameDaughters)
  ,_lctTensor(LeviCivitaTensor())
  ,_metricTensor(MetricTensor())
  ,_polMother(PolVector(motherJPC->J))
   ,_polDaughter1(PolVector(daughter1JPC->J))
   ,_polDaughter2(PolVector(daughter2JPC->J))
  ,_qCMSReal(0.)
  ,_mother4Vec(1., 0.,0.,0.)
{
  _theTFile=new TFile("DecAngularDist.root","recreate");
  _allLSs=GetValidLS(_motherJPC->J, _motherJPC->P, _daughter1JPC->J,_daughter1JPC->P, _daughter2JPC->J, _daughter2JPC->P);


  //check C parity
  bool allowedC=true;
  if (!sameDaughters){
    if (_motherJPC->C != (daughter1JPC->C * daughter2JPC->C)){
      Alert << "no C-parity conservation for this decay" << endmsg;
      exit(1);
    }
  }

  std::vector<LS>::const_iterator itLS;

  for (itLS=_allLSs.begin(); itLS!=_allLSs.end(); ++itLS){

    if (sameDaughters){
      bool bosons=true;
      if (int(2*daughter1JPC->J)%2 !=0) bosons=false;
      
      Spin theL=itLS->L;
      Spin theS=itLS->S;
      int cparity=0;
      
      if (bosons) cparity=int(pow(-1,theL));
      else cparity=int(pow(-1,(theL+theS))); //daughters are fermions	
      
      if(_motherJPC->C != cparity ) allowedC=false;
    }
    
    if (allowedC) _allowedCParityLSs.push_back(*itLS);
    else _forbiddenCParityLSs.push_back(*itLS);
  }

 for (Spin lamMother=-_motherJPC->J; lamMother<=_motherJPC->J; ++lamMother){
   std::stringstream lamMotherStrStr;
   lamMotherStrStr << lamMother;

   std::string histThetaTensorAmpJMName="ampJMThetaTensor"+_motherJPC->name()+"_"+lamMotherStrStr.str();

   for (Spin lamDaughter1=-_daughter1JPC->J; lamDaughter1<=_daughter1JPC->J; ++lamDaughter1){
     std::stringstream M1StrStr;
     M1StrStr << lamDaughter1;
     for (Spin lamDaughter2=-_daughter2JPC->J; lamDaughter2<=_daughter2JPC->J; ++lamDaughter2){
       std::stringstream M2StrStr;
       M2StrStr << lamDaughter2;

       std::string histThetaTensorJMM1M1Name="TensorJMM1M2"+_motherJPC->name()+"_"+lamMotherStrStr.str()+"_M1_"+M1StrStr.str()+"_M2_"+M2StrStr.str();
       TH1F* currentTensorJMM1M2Hist=new TH1F(histThetaTensorJMM1M1Name.c_str(), histThetaTensorJMM1M1Name.c_str(), 202., -1.01, 1.01);
       _tensorJMM1M2AmpThetaHistMap[lamMother][lamDaughter1][lamDaughter2]=currentTensorJMM1M2Hist;

       std::string histThetaCanoJMM1M1Name="CanoJMM1M2"+_motherJPC->name()+"_"+lamMotherStrStr.str()+"_M1_"+M1StrStr.str()+"_M2_"+M2StrStr.str();
       TH1F* currentCanoJMM1M2Hist=new TH1F(histThetaCanoJMM1M1Name.c_str(), histThetaCanoJMM1M1Name.c_str(), 202., -1.01, 1.01);
       _canoJMM1M2AmpThetaHistMap[lamMother][lamDaughter1][lamDaughter2]=currentCanoJMM1M2Hist;

       for (itLS=_allowedCParityLSs.begin(); itLS!=_allowedCParityLSs.end(); ++itLS){
	 std::stringstream lStrStr;
	 lStrStr << (*itLS).L;
	 std::stringstream sStrStr;
	 sStrStr << (*itLS).S;

	 std::string histThetaTensorJLSMM1M1Name="Tensor"+_motherJPC->name()+"_L"+lStrStr.str()+"_S"+sStrStr.str()+"_M"+lamMotherStrStr.str()+"_M1_"+M1StrStr.str()+"_M2_"+M2StrStr.str();
	 TH1F* currentTensorJLSMM1M2Hist=new TH1F(histThetaTensorJLSMM1M1Name.c_str(), histThetaTensorJLSMM1M1Name.c_str(), 202., -1.01, 1.01);
	 _tensorJLSMM1M2AmpThetaHistMap[(*itLS).L][(*itLS).S][lamMother][lamDaughter1][lamDaughter2]=currentTensorJLSMM1M2Hist;

	 std::string histThetaCanoJLSMM1M1Name="Cano"+_motherJPC->name()+"_L"+lStrStr.str()+"_S"+sStrStr.str()+"_M"+lamMotherStrStr.str()+"_M1_"+M1StrStr.str()+"_M2_"+M2StrStr.str();
	 TH1F* currentCanoJLSMM1M2Hist=new TH1F(histThetaCanoJLSMM1M1Name.c_str(), histThetaCanoJLSMM1M1Name.c_str(), 202., -1.01, 1.01);
	 _canoJLSMM1M2AmpThetaHistMap[(*itLS).L][(*itLS).S][lamMother][lamDaughter1][lamDaughter2]=currentCanoJLSMM1M2Hist;
	 }
     }
   }
 }
 
 _heliJAmpAllThetaHist = new TH1F("ampJAllThetaHeli", "ampJAllThetaHeli", 202., -1.01, 1.01);
 _heliLSJAmpAllThetaHist = new TH1F("ampJAllThetaHeliLS", "ampJAllThetaHeliLS", 202., -1.01, 1.01);
 _canoJAmpAllThetaHist = new TH1F("ampJAllThetaCano", "ampJAllThetaCano", 202., -1.01, 1.01);
 _tensorJAmpAllThetaHist = new TH1F("ampJAllThetaTensor", "ampJAllThetaTensor", 202., -1.01, 1.01);
 _tensorJMM1M2AmpAllThetaHist = new TH1F("ampJMM1M2AllThetaTensor", "ampJMM1M2AllThetaTensor", 202., -1.01, 1.01);
 _canoJMM1M2AmpAllThetaHist = new TH1F("ampJMM1M2AllThetaCano", "ampJMM1M2AllThetaCano", 202., -1.01, 1.01);

  //fill heli
  _polMother.SetP4(_mother4Vec, _mother4Vec.Mass());
  complex<double> qCMS=breakupMomQ(1., .135, .135);
  _qCMSReal=qCMS.real();


  for (double cosT=-1.; cosT<=1.; cosT+=2./200.){
    fillCanoAmps(cosT);
    fillheliAmps(cosT);
    fillTensorAmps(cosT);
  }
}

DecAngularDist::~DecAngularDist()
{
  _theTFile->Write();
  _theTFile->Close();
}

void DecAngularDist::calcChi12(Tensor<complex<double> >& s12SpinProjector, Tensor<complex<double> >& epsilonDaughter1Project, Tensor<complex<double> >& epsilonDaughter2Project, Vector4<double>& mother_4Vec, Tensor<complex<double> >& result){
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
 // DebugMsg << "leviPssTensor: " << leviPssTensor << endmsg;

  Tensor<complex<double> > chiPart;
  if(add_lctForChi){
    chiPart=leviPssTensor*conj(epsilonDaughter1Project);
    //    chiPart.Permute(1,chiPart.Rank());
    chiPart.Permute(0,chiPart.Rank()-1);
    chiPart=chiPart*conj(epsilonDaughter2Project);
  }
  else{
    DebugMsg << "conj(epsilonDaughter1Project):\t" << conj(epsilonDaughter1Project) << endmsg;
    DebugMsg << "conj(epsilonDaughter2Project):\t" << conj(epsilonDaughter2Project) << endmsg;
    if(spin12==0) chiPart=conj(epsilonDaughter1Project)|conj(epsilonDaughter2Project);
    else chiPart=conj(epsilonDaughter1Project)%conj(epsilonDaughter2Project);
  }
  //  DebugMsg << "chiPart:\t" << chiPart << endmsg;

  int noContractionsChi12=0.5*(spin12+spin1+spin2);
  int noRemainingContr=noContractionsChi12;
  if (add_lctForChi){
    noContractionsChi12=0.5*(spin12+4+1+spin1+spin2);
    noRemainingContr=noContractionsChi12-3;    
  }

  result=s12SpinProjector.Contract(chiPart, noRemainingContr);
  //  DebugMsg << "no of contractions chi12:\t" << noContractionsChi12 << endmsg;
  int rankCh12=2*spin12+spin1+spin2-2*noContractionsChi12;
  if(add_lctForChi) rankCh12=2*spin12+spin1+spin2+5-2*noContractionsChi12; 
  // DebugMsg << "chi12:\t" << "should have rank " << rankCh12 << " and (cross check): " << spin12 <<"\trank" << result.Rank() << endmsg;
  // if(result.Rank()<3){
  //   DebugMsg << result << endmsg;
  // }
}


void DecAngularDist::calcLSpart(OrbitalTensor& orbTensor, Tensor<complex<double> >& chi12, Spin spin12, Vector4<double>& mother_4Vec, Tensor<complex<double> >& result){
  int orbitalL=orbTensor.Rank();
  // if(orbTensor.Rank()<3){
  //   DebugMsg << "orbTensor: " << orbTensor << endmsg;    
  // }


  int s12LJ=spin12+orbitalL+_motherJPC->J; 
  Tensor<complex<double> > leviPlsTensor;
  bool add_lct=true;
  int noOfContractions=0;
  if( s12LJ%2 ==0 ){ //even
    add_lct=false;
    leviPlsTensor=complex<double>(1.,0.);
    noOfContractions=(orbitalL+spin12-_motherJPC->J)/2;
  }
  else{
    leviPlsTensor=_lctTensor*mother_4Vec;
    noOfContractions=(4+1+orbitalL+spin12-_motherJPC->J)/2;
  }
  // DebugMsg << "leviPlsTensor: " << leviPlsTensor << endmsg;
  // DebugMsg << "no of contr.: " << noOfContractions << endmsg;

  Tensor<complex<double> > lsPartTensor;
  if(add_lct){
    lsPartTensor=leviPlsTensor*orbTensor;
    // lsPartTensor.Permute(1,lsPartTensor.Rank());
    lsPartTensor.Permute(0,lsPartTensor.Rank()-1);
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

  //    DebugMsg << "LS part: should have rank " << _motherJPC->J << "\t rank" << result.Rank() << endmsg;
    // if(result.Rank()<3){
    //   DebugMsg << result << endmsg;
    // }
}


void  DecAngularDist::fillTensorAmps(double cosT){
	double theta = acos(cosT);

	TLorentzVector daughter1LV( 0., 0., _qCMSReal, sqrt(.135*.135+_qCMSReal*_qCMSReal));
	daughter1LV.SetTheta(theta);
	Vector4<double> daughter1_4Vec(daughter1LV.E(), daughter1LV.Px(), daughter1LV.Py(), daughter1LV.Pz());
	Vector4<double> daughter2_4Vec(daughter1LV.E(), -daughter1LV.Px(), -daughter1LV.Py(), -daughter1LV.Pz());

	std::cout << "daughter1_4Vec.Phi(): " << daughter1_4Vec.Phi() << std::endl;
	std::cout << "daughter2_4Vec.Phi(): " << daughter2_4Vec.Phi() << std::endl;

	_polDaughter1.SetP4(daughter1_4Vec, daughter1_4Vec.Mass());
	_polDaughter2.SetP4(daughter2_4Vec, daughter2_4Vec.Mass());

	for (Spin lamMother=-_motherJPC->J; lamMother<=_motherJPC->J; ++lamMother){
	  complex<double> sumAmp(0.,0.);
	  double sumInCoherentAmp=0.;
	  Tensor<complex<double> > epsilonMotherProject = _polMother(lamMother);

	  std::map< Spin, std::map< Spin, complex<double> > > ampM1M2Map;
	  for (Spin lamDaughter1=-_daughter1JPC->J; lamDaughter1<=_daughter1JPC->J; ++lamDaughter1){
	    for (Spin lamDaughter2=-_daughter2JPC->J; lamDaughter2<=_daughter2JPC->J; ++lamDaughter2){
	      ampM1M2Map[lamDaughter1][lamDaughter2]=complex<double>(0.,0.);
	    }
	  }


	  std::vector<LS>::const_iterator itLS;	  
	  for (itLS=_allowedCParityLSs.begin(); itLS!=_allowedCParityLSs.end(); ++itLS){

	    complex<double> currentTensorLSAmp(0.,0.);
	    OrbitalTensor orbTensor((*itLS).L);
	    orbTensor.SetP4(daughter1_4Vec, daughter2_4Vec);
	    
	    //calculate chi
	    Tensor<complex<double> > s12SpinProjector;
	    PolVector part12PolVec((*itLS).S);
	    part12PolVec.SetP4(_mother4Vec, _mother4Vec.Mass());

	    s12SpinProjector=part12PolVec.Projector();

	    for (Spin lamDaughter1=-_daughter1JPC->J; lamDaughter1<=_daughter1JPC->J; ++lamDaughter1){
	      Tensor<complex<double> > epsilonDaughter1Project = _polDaughter1(lamDaughter1);

	      for (Spin lamDaughter2=-_daughter2JPC->J; lamDaughter2<=_daughter2JPC->J; ++lamDaughter2){
		//		Spin lamDaughterDiff=lamDaughter1-lamDaughter2;
		//		if( fabs(lamDaughterDiff)>_motherJPC->J || fabs(lamDaughterDiff)>(*itLS).S) continue;
	
		Tensor<complex<double> > epsilonDaughter2Project = _polDaughter2(lamDaughter2);
		Tensor<complex<double> > chi12;
		calcChi12(s12SpinProjector, epsilonDaughter1Project, epsilonDaughter2Project, _mother4Vec, chi12);
		// calcChi12(s12SpinProjector, epsilonDaughter1Project, epsilonDaughter2Project,mother4VecBoost, chi12);		
		
		if (chi12.Rank()<2) std::cout << "chi12: " << chi12 << std::endl;
		
		//calculate ls part;
		Tensor<complex<double> > lsPartTensorNew;
		calcLSpart(orbTensor, chi12, (*itLS).S, _mother4Vec, lsPartTensorNew);
		//		calcLSpart(orbTensor, chi12, (*itLS).S, mother4VecBoost, lsPartTensorNew);		
		
		std::cout << "lsPartTensorNew.Rank(): " << lsPartTensorNew.Rank() << std::endl;		
		if (lsPartTensorNew.Rank()<3) std::cout << "lsPartTensorNew: " << lsPartTensorNew << std::endl;
		
		Tensor<complex<double> > finalTmpAmp=epsilonMotherProject|lsPartTensorNew;
		currentTensorLSAmp+= finalTmpAmp(0);
		sumAmp+=finalTmpAmp(0);
		sumInCoherentAmp+=norm(finalTmpAmp(0));
		ampM1M2Map.at(lamDaughter1).at(lamDaughter2)+=finalTmpAmp(0);
		_tensorJLSMM1M2AmpThetaHistMap.at((*itLS).L).at((*itLS).S).at(lamMother).at(lamDaughter1).at(lamDaughter2)->Fill(cosT, finalTmpAmp(0).real());
		//_tensorJLSMM1M2AmpThetaHistMap.at((*itLS).L).at((*itLS).S).at(lamMother).at(lamDaughter1).at(lamDaughter2)->Fill(cosT, norm(finalTmpAmp(0)));
	      }
	    }
	  }

	  for (Spin lamDaughter1=-_daughter1JPC->J; lamDaughter1<=_daughter1JPC->J; ++lamDaughter1){
	    for (Spin lamDaughter2=-_daughter2JPC->J; lamDaughter2<=_daughter2JPC->J; ++lamDaughter2){
	     _tensorJMM1M2AmpThetaHistMap.at(lamMother).at(lamDaughter1).at(lamDaughter2)->Fill(cosT, norm(ampM1M2Map.at(lamDaughter1).at(lamDaughter2)));
	     // _tensorJMM1M2AmpThetaHistMap.at(lamMother).at(lamDaughter1).at(lamDaughter2)->Fill(cosT, ampM1M2Map.at(lamDaughter1).at(lamDaughter2).real());
	      _tensorJMM1M2AmpAllThetaHist->Fill(cosT, norm(ampM1M2Map.at(lamDaughter1).at(lamDaughter2)));
	    }
	  }
	  _tensorJAmpAllThetaHist->Fill(cosT, sumInCoherentAmp);
	} 
	
	
}


void DecAngularDist::fillheliAmps(double cosT){
  double theta=acos(cosT);
  Spin sMax=0;
  std::vector<LS>::const_iterator itLS;  
  for (itLS=_allowedCParityLSs.begin(); itLS!=_allowedCParityLSs.end(); ++itLS){
    if(sMax<(*itLS).S) sMax=(*itLS).S; 
  }

  for (Spin lamMother=-_motherJPC->J; lamMother<=_motherJPC->J; ++lamMother){
    for (Spin lamDaughter1=-_daughter1JPC->J; lamDaughter1<=_daughter1JPC->J; ++lamDaughter1){
      for (Spin lamDaughter2=-_daughter2JPC->J; lamDaughter2<=_daughter2JPC->J; ++lamDaughter2){
	complex<double> ampSum(0.,0.);
	complex<double> ampLSSum(0.,0.);
	Spin lamDaughterDiff=lamDaughter1-lamDaughter2;
	if(fabs(lamDaughterDiff)>_motherJPC->J || fabs(lamDaughterDiff)>sMax) continue;
	ampSum=  conj(Wigner_D(0., theta,0, _motherJPC->J, lamMother, lamDaughterDiff));
	
	for (itLS=_allowedCParityLSs.begin(); itLS!=_allowedCParityLSs.end(); ++itLS){
	  if(fabs(lamDaughterDiff)>_motherJPC->J || fabs(lamDaughterDiff)>(*itLS).S) continue;
	  ampLSSum+=sqrt(2.*(*itLS).L+1)
	    *Clebsch((*itLS).L, 0, (*itLS).S, lamDaughterDiff, _motherJPC->J, lamDaughterDiff)
	    *Clebsch(_daughter1JPC->J, lamDaughter1, _daughter2JPC->J, -lamDaughter2, (*itLS).S, lamDaughterDiff)
	    *conj(Wigner_D(0., theta,0, _motherJPC->J, lamMother, lamDaughterDiff));
	}

	_heliJAmpAllThetaHist->Fill(cosT, norm(ampSum));
	
	_heliLSJAmpAllThetaHist->Fill(cosT, norm(ampLSSum));
	
      }
    }
    
  }
}

void  DecAngularDist::fillCanoAmps(double cosT){
  double theta = acos(cosT);
  
  for (Spin lamMother=-_motherJPC->J; lamMother<=_motherJPC->J; ++lamMother){
    complex<double> currentJMamp(0.,0.);
    complex<double> currentJMLSamp(0.,0.);
    
    std::map< Spin, std::map< Spin, complex<double> > > ampM1M2Map;
    for (Spin lamDaughter1=-_daughter1JPC->J; lamDaughter1<=_daughter1JPC->J; ++lamDaughter1){
      for (Spin lamDaughter2=-_daughter2JPC->J; lamDaughter2<=_daughter2JPC->J; ++lamDaughter2){
	ampM1M2Map[lamDaughter1][lamDaughter2]=complex<double>(0.,0.);
      }
    }
    
    std::vector<LS>::const_iterator itLS;         
    for (itLS=_allowedCParityLSs.begin(); itLS!=_allowedCParityLSs.end(); ++itLS){
      complex <double> currentLSAmp(0., 0.);
      complex <double> currentLSHeliAmp(0., 0.);
      
      for (Spin mDaughter1=-_daughter1JPC->J; mDaughter1<=_daughter1JPC->J; ++mDaughter1){
	for (Spin mDaughter2=-_daughter2JPC->J; mDaughter2<=_daughter2JPC->J; ++mDaughter2){
	  
	  Spin mL=lamMother-mDaughter1-mDaughter2;
	  if( fabs(mL)>(*itLS).L ) continue;
	  
	  Spin mS=mDaughter1+mDaughter2;
	  if( fabs(mS)>(*itLS).S ) continue;
	  
	  complex<double> sphereHarmAmp=SphericalHarmonic::SphericalHarmonicFunc((*itLS).L, mL, theta, 0.);
	  double cgLS=Clebsch((*itLS).L, mL, (*itLS).S, mS, _motherJPC->J, lamMother);
	  
	  double cgJJ=Clebsch(_daughter1JPC->J, mDaughter1, _daughter2JPC->J, mDaughter2, (*itLS).S, mS);
	  
	  complex<double> currentMM1M2Amp=cgLS*cgJJ*sphereHarmAmp;
	  ampM1M2Map.at(mDaughter1).at(mDaughter2)+=currentMM1M2Amp;
	  _canoJLSMM1M2AmpThetaHistMap.at((*itLS).L).at((*itLS).S).at(lamMother).at(mDaughter1).at(mDaughter2)->Fill(cosT, norm(currentMM1M2Amp));
	}
      }
    }
    
    for (Spin lamDaughter1=-_daughter1JPC->J; lamDaughter1<=_daughter1JPC->J; ++lamDaughter1){
      for (Spin lamDaughter2=-_daughter2JPC->J; lamDaughter2<=_daughter2JPC->J; ++lamDaughter2){
	//	_canoJMM1M2AmpThetaHistMap.at(lamMother).at(lamDaughter1).at(lamDaughter2)->Fill(cosT, norm(ampM1M2Map.at(lamDaughter1).at(lamDaughter2)));
      	_canoJMM1M2AmpThetaHistMap.at(lamMother).at(lamDaughter1).at(lamDaughter2)->Fill(cosT, ampM1M2Map.at(lamDaughter1).at(lamDaughter2).real());
	_canoJMM1M2AmpAllThetaHist->Fill(cosT, norm(ampM1M2Map.at(lamDaughter1).at(lamDaughter2)));

      }
    }

  }
}




