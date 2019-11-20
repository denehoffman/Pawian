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

#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/AbsStates.hh"
#include "pbarpUtils/pbarpStatesLS.hh"
#include "PwaUtils/DataUtils.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

int main(int __argc,char *__argv[]){

  if( __argc>1 && ( strcmp( __argv[1], "-help" ) == 0
		    || strcmp( __argv[1], "--help" ) == 0 ) ){

    std::cout << "test application which calculates all initial pbar p amplitudes in the tensor formalism\n"
	      << "with -msg <errorLogMode>: choose the mode for the error logger \n"
	      << "with -lmax <value>: choose the maximum Spin J \n"
	      << std::endl;
    return 0;
  }

  int optind=1;
  std::string msgModeStr="default";
  std::string lmaxStr="0";
  // decode arguments

  while ((optind < (__argc-1) ) && (__argv[optind][0]=='-')) {
    bool found=false;
    std::string sw = __argv[optind];
    if (sw=="-msg"){
      optind++;
      msgModeStr = __argv[optind];
      found=true;
    }
    if(sw=="-lmax"){
      optind++;
      lmaxStr = __argv[optind];
      found=true;
    }
    if (!found){
      WarningMsg << "Unknown switch: " << __argv[optind] << endmsg;
      optind++;
    }
    
   while ( (optind < __argc ) && __argv[optind][0]!='-' ) optind++;
  }

  if(msgModeStr == "debugging") ErrLogger::instance().setThreshold(logging::log_level::DEBUG);
  else if(msgModeStr == "trace") ErrLogger::instance().setThreshold(logging::log_level::INFO);
  else if(msgModeStr == "routine") ErrLogger::instance().setThreshold(logging::log_level::INFO);
  else if(msgModeStr == "warning") ErrLogger::instance().setThreshold(logging::log_level::WARN);
  else if(msgModeStr == "error")   ErrLogger::instance().setThreshold(logging::log_level::ERROR);
  else {
    ErrLogger::instance().setThreshold(logging::log_level::INFO);
    WarningMsg << "ErrorLogger not (properly) set -> Use mode 'WARN' " << endmsg;  
  }

 std::stringstream lmaxStrStr(lmaxStr);
 int lmax=2;
 lmaxStrStr >> lmax;

  InfoMsg << "lmax: " << lmax << endmsg;

  pbarpStatesLS thepbarpState(lmax);

  thepbarpState.print(std::cout);

  double pbarMass=0.93827;
  double pMass=0.93827;
  double pbarMom=0.9;

  Vector4<double> pbar4Vec(sqrt(pbarMass*pbarMass+pbarMom*pbarMom), 0.,0., pbarMom);
  Vector4<double> p4Vec(pMass, 0.,0.,0.);
  Vector4<double> allVec=pbar4Vec+p4Vec;

  std::cout << "\n\npbar4Vec: " <<  pbar4Vec << std::endl;
  std::cout << "p4Vec: " <<  p4Vec << std::endl;
  std::cout << "allVec: " <<  allVec << std::endl;

  Vector4<double> pbarCM4Vec=pbar4Vec;
  Vector4<double> pCM4Vec=p4Vec;
  Vector4<double> allCM4Vec=allVec;

  pbarCM4Vec.Boost(allVec);
  pCM4Vec.Boost(allVec);
  allCM4Vec.Boost(allVec);

  std::cout << "\npbarCM4Vec: " <<  pbarCM4Vec << std::endl;
  std::cout << "pCM4Vec: " <<  pCM4Vec << std::endl;
  std::cout << "allCM4Vec: " <<  allCM4Vec << std::endl;  

  DiracGamma5 gamma5;
   DiracGamma gamma;

  DiracAntiSpinor pbarAntiSpinor; // pbar anti spinor
  DiracSpinor pbarSpinor; // pbar spinor
  DiracSpinor pSpinor; // p spinor

  pbarAntiSpinor.SetP4(pbar4Vec, pbarMass);
  pbarSpinor.SetP4(pbar4Vec, pbarMass);
  pSpinor.SetP4(p4Vec, pMass);

  // std::cout << "\npbarAntiSpinor(1/2):\n" <<  pbarAntiSpinor(1./2) << std::endl;
  // std::cout << "\nBar(pbarAntiSpinor(1/2)):\n" <<  Bar(pbarAntiSpinor(1./2)) << std::endl;
  // std::cout << "\npbarAntiSpinor(-1/2):\n" <<  pbarAntiSpinor(-1./2) << std::endl;
  // std::cout << "\nBar(pbarAntiSpinor(-1/2)):\n" <<  Bar(pbarAntiSpinor(-1./2)) << std::endl;
  // std::cout << "pbarSpinor(1/2):\n" <<  pbarSpinor(1./2) << std::endl;
  // std::cout << "Bar(pbarSpinor(1/2)):\n" <<  Bar(pbarSpinor(1./2)) << std::endl;
  // std::cout << "pbarSpinor(-1/2):\n" <<  pbarSpinor(-1./2) << std::endl;
  // std::cout << "Bar(pbarSpinor(-1/2)):\n" <<  Bar(pbarSpinor(-1./2)) << std::endl;
  // std::cout << "pSpinor(1/2):\n" <<  pSpinor(1./2) << std::endl; 
  // std::cout << "pSpinor(-1/2):\n" <<  pSpinor(-1./2) << std::endl; 
  

  DiracAntiSpinor pbarAntiSpinorCM; // pbar spinor
  DiracSpinor pbarSpinorCM; // pbar spinor
  DiracSpinor pSpinorCM; // p spinor

  pbarAntiSpinorCM.SetP4(pbarCM4Vec, pbarMass);
  pbarSpinorCM.SetP4(pbarCM4Vec, pbarMass);
  pSpinorCM.SetP4(pCM4Vec, pMass);

  // std::cout << "\npbarAntiSpinorCM(1/2):\n" <<  pbarAntiSpinorCM(1./2) << std::endl;
  // std::cout << "\npbarAntiSpinorCM(-1/2):\n" <<  pbarAntiSpinorCM(-1./2) << std::endl;
  // std::cout << "pbarSpinorCM(1/2):\n" <<  pbarSpinorCM(1./2) << std::endl;
  // std::cout << "pbarSpinorCM(-1/2):\n" <<  pbarSpinorCM(-1./2) << std::endl;
  // std::cout << "pSpinorCM(1/2):\n" <<  pSpinorCM(1./2) << std::endl; 
  // std::cout << "pSpinorCM(-1/2):\n" <<  pSpinorCM(-1./2) << std::endl;


  double psi0AtRestNom=(pbarCM4Vec.E()+pbarCM4Vec.M())*(pCM4Vec.E()+pCM4Vec.M())+pCM4Vec.P()*pCM4Vec.P();
  double psi0AtRestDemom=sqrt(4.*pbarCM4Vec.M()*pCM4Vec.M()*(pbarCM4Vec.E()+pbarCM4Vec.M())*(pCM4Vec.E()+pCM4Vec.M()));
  double psi0AtRest=psi0AtRestNom/psi0AtRestDemom;
  std::cout << "\nPsi0AtRest:\n" <<  psi0AtRest << std::endl;


  double psi1Fac1=sqrt((pbarCM4Vec.E()+pbarCM4Vec.M())*(pCM4Vec.E()+pCM4Vec.M())/(4.*pbarCM4Vec.M()*pCM4Vec.M()));
  double psi1Fac2=1.+pCM4Vec.P()*pCM4Vec.P()/((pbarCM4Vec.E()+pbarCM4Vec.M())*(pCM4Vec.E()+pCM4Vec.M()));

  Matrix<complex<double> > chiP12(2,1);
  chiP12(0,0) = 1.;  
  chiP12(1,0) = 0.;     
  Matrix<complex<double> > chiPm12(2,1);
  chiPm12(0,0) = 0.;  
  chiPm12(1,0) = 1.;

  Matrix<complex<double> > chiPbar12Dagger(2,1); 
  chiPbar12Dagger=chiP12.Adjoint();

  Matrix<complex<double> > chiPbarm12Dagger(2,1); 
  chiPbarm12Dagger=chiPm12.Adjoint();

  PauliSigma pauliSigma;

  Matrix<complex<double> > psi1_12_m12_0=(chiPbar12Dagger*pauliSigma[0]*chiPm12);
  std::cout << " psi1_12_m12_0: " <<  psi1_12_m12_0 << std::endl;

  Matrix<complex<double> > psi1_12_m12_1=(chiPbar12Dagger*pauliSigma[1]*chiPm12);
  std::cout << " psi1_12_m12_1: " <<  psi1_12_m12_1 << std::endl;
  
  Matrix<complex<double> > psi1_12_12_0=psi1Fac1*psi1Fac2*(chiPbar12Dagger*pauliSigma[0]*chiP12);
  std::cout << " psi1_12_12_0: " <<  psi1_12_12_0 << std::endl;

  Matrix<complex<double> > psi1_12_12_1=psi1Fac1*psi1Fac2*(chiPbar12Dagger*pauliSigma[1]*chiP12);
  std::cout << " psi1_12_12_1: " <<  psi1_12_12_1 << std::endl;

  Matrix<complex<double> > psi1_m12_m12_0=psi1Fac1*psi1Fac2*(chiPbarm12Dagger*pauliSigma[0]*chiPm12);
  std::cout << " psi1_m12_m12_0: " <<  psi1_m12_m12_0 << std::endl;

  Matrix<complex<double> > psi1_m12_m12_1=psi1Fac1*psi1Fac2*(chiPbarm12Dagger*pauliSigma[1]*chiPm12);
  std::cout << " psi1_m12_m12_1: " <<  psi1_m12_m12_1 << std::endl;

  std::cout << "\nBar(pSpinorCM(1./2))*gamma5*pbarAntiSpinorCM(-1./2):\n" <<  Bar(pSpinorCM(1./2))*gamma5*pbarAntiSpinorCM(-1./2) << std::endl;
  std::cout << "\nBar(pSpinor(1./2))*gamma5*pbarAntiSpinor(-1./2):\n" <<  Bar(pSpinor(1./2))*gamma5*pbarAntiSpinor(-1./2) << std::endl;

  Tensor<complex<double> > Psi0_m12_12=Bar(pbarAntiSpinor(1./2))*gamma5*pSpinor(-1./2);
  Tensor<complex<double> > Psi0_m12_12CM=Bar(pbarAntiSpinorCM(-1./2))*gamma5*pSpinorCM(1./2);

  std::cout << "\nPsi0_m12_12:\n" <<  Psi0_m12_12 << std::endl;
  std::cout << "\nPsi0_m12_12CM:\n" <<  Psi0_m12_12CM << std::endl;

  Tensor<complex<double> > Psi0_12_m12=Bar(pbarAntiSpinor(1./2))*gamma5*pSpinor(-1./2);
  Tensor<complex<double> > Psi0_12_m12CM=Bar(pbarAntiSpinorCM(1./2))*gamma5*pSpinorCM(-1./2);

  std::cout << "\nPsi0_12_m12:\n" <<  Psi0_12_m12 << std::endl;
  std::cout << "\nPsi0_12_m12CM:\n" <<  Psi0_12_m12CM << std::endl;


  // Tensor<complex<double> > Psi0_12_12=Bar(pSpinor(1./2))*gamma5*pbarAntiSpinor(1./2);
  // Tensor<complex<double> > Psi0_12_12CM=Bar(pSpinorCM(1./2))*gamma5*pbarAntiSpinorCM(1./2);

  // std::cout << "\nPsi0_12_12:\n" <<  Psi0_12_12 << std::endl;
  // std::cout << "\nPsi0_12_12CM:\n" <<  Psi0_12_12CM << std::endl;

  // Tensor<complex<double> > Psi0_m12_m12=Bar(pSpinor(-1./2))*gamma5*pbarAntiSpinor(-1./2);
  // Tensor<complex<double> > Psi0_m12_m12CM=Bar(pSpinorCM(-1./2))*gamma5*pbarAntiSpinorCM(-1./2);

  // std::cout << "\nPsi0_m12_m12:\n" <<  Psi0_m12_m12 << std::endl;
  // std::cout << "\nPsi0_m12_m12CM:\n" <<  Psi0_m12_m12CM << std::endl;

  Tensor<complex<double> > cancelTerm=(pbar4Vec-p4Vec)/(pbarMass+pMass+allVec.M());
  Tensor<complex<double> > cancelTermCM=(pbarCM4Vec-pCM4Vec)/(pbarMass+pMass+allCM4Vec.M());

  Tensor<complex<double> > Psi1_12_12=Bar(pbarAntiSpinor(1./2))*gamma*pSpinor(1./2)-Bar(pbarAntiSpinor(1./2))*cancelTerm*pSpinor(1./2);
  Tensor<complex<double> > Psi1_12_12CM=Bar(pbarAntiSpinorCM(1./2))*gamma*pSpinorCM(1./2)-Bar(pbarAntiSpinorCM(1./2))*cancelTermCM*pSpinor(1./2);

  std::cout << "\nBar(pSpinor(1./2))*gamma*pAntiSpinor(1./2):" << Bar(pSpinor(1./2))*gamma*pbarAntiSpinor(1./2) << std::endl;
  std::cout << "\nBar(pSpinor(-1./2))*gamma*pAntiSpinor(-1./2):" << Bar(pSpinor(-1./2))*gamma*pbarAntiSpinor(-1./2) << std::endl;
  std::cout << "\nBar(pSpinor(1./2))*gamma*pAntiSpinor(-1./2):" << Bar(pSpinor(1./2))*gamma*pbarAntiSpinor(-1./2) << std::endl;

  std::cout << "\nBar(pSpinor(1./2))*cancelTerm*pbarAntiSpinor(1./2):" << Bar(pSpinor(1./2))*cancelTerm*pbarAntiSpinor(1./2) << std::endl;
  std::cout << "\nBar(pbarAntiSpinor(1./2))*cancelTerm*pSpinor(1./2):" << Bar(pbarAntiSpinor(1./2))*cancelTerm*pSpinor(1./2) << std::endl;

  std::cout << "\nBar(pbarSpinor(1./2))*cancelTerm*pSpinor(1./2):" << Bar(pbarSpinor(1./2))*cancelTerm*pSpinor(1./2) << std::endl;
  std::cout << "\nBar(pbarAntiSpinor(-1./2))*cancelTerm*pSpinor(1./2):" << Bar(pbarAntiSpinor(-1./2))*cancelTerm*pSpinor(1./2) << std::endl;
  std::cout << "\nBar(pbarAntiSpinor(1./2))*cancelTerm*pSpinor(-1./2):" << Bar(pbarAntiSpinor(-1./2))*cancelTerm*pSpinor(1./2) << std::endl;
  //  std::cout << "\npbarSpinor(1./2)*cancelTerm*pSpinor(1./2):" << pbarSpinor(1./2)*cancelTerm*pSpinor(1./2) << std::endl;

  std::cout << "\nPsi1_12_12:\n" <<  Psi1_12_12 << std::endl;
  std::cout << "\nPsi1_12_12CM:\n" <<  Psi1_12_12CM << std::endl;

  Tensor<complex<double> > Psi1_m12_12=Bar(pbarAntiSpinor(-1./2))*gamma*pSpinor(1./2)-Bar(pbarAntiSpinor(-1./2))*cancelTerm*pSpinor(1./2);
  Tensor<complex<double> > Psi1_m12_12CM=Bar(pbarAntiSpinorCM(-1./2))*gamma*pSpinorCM(1./2)-Bar(pbarAntiSpinorCM(-1./2))*cancelTermCM*pSpinorCM(1./2);

  std::cout << "\nPsi1_m12_12:\n" <<  Psi1_m12_12 << std::endl;
  std::cout << "\nPsi1_m12_12CM:\n" <<  Psi1_m12_12CM << std::endl;

  Tensor<complex<double> > Psi1_12_m12=Bar(pbarAntiSpinor(1./2))*gamma*pSpinor(-1./2)-Bar(pbarAntiSpinor(1./2))*cancelTerm*pSpinor(-1./2);
  Tensor<complex<double> > Psi1_12_m12CM=Bar(pbarAntiSpinorCM(1./2))*gamma*pSpinorCM(-1./2)-Bar(pbarAntiSpinorCM(1./2))*cancelTermCM*pSpinorCM(-1./2);

  std::cout << "\nPsi1_12_m12:\n" <<  Psi1_12_m12 << std::endl;
  std::cout << "\nPsi1_12_m12CM:\n" <<  Psi1_12_m12CM << std::endl;

  Tensor<complex<double> > Psi1_m12_m12=Bar(pbarAntiSpinor(-1./2))*gamma*pSpinor(-1./2)-Bar(pbarAntiSpinor(-1./2))*cancelTerm*pSpinor(-1./2);
  Tensor<complex<double> > Psi1_m12_m12CM=Bar(pbarAntiSpinorCM(-1./2))*gamma*pSpinorCM(-1./2)-Bar(pbarAntiSpinorCM(-1./2))*cancelTermCM*pSpinorCM(-1./2);

  std::cout << "\nPsi1_m12_m12:\n" <<  Psi1_m12_m12 << std::endl;
  std::cout << "\nPsi1_m12_m12CM:\n" <<  Psi1_m12_m12CM << std::endl;

  std::cout <<"\n\n\n" << "*** singlet states ******" <<"\n\n\n" << std::endl;
  std::vector< std::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< std::shared_ptr<const JPCLS> > jpclsStatesSiglet=thepbarpState.singlet_JPCLS_States();

  for(itJPCLS=jpclsStatesSiglet.begin(); itJPCLS!=jpclsStatesSiglet.end(); ++itJPCLS){
    OrbitalTensor currentTensor((*itJPCLS)->L);
    currentTensor.SetP4(pbar4Vec,p4Vec);

    PolVector currentPolVec((*itJPCLS)->J);
    currentPolVec.SetP4(allVec, allVec.M());

    Tensor<complex<double> > result= currentPolVec(0).Conjugate() | currentTensor;
    result*= Psi0_m12_12;
     
    std::cout <<"\n" << (*itJPCLS)->name() << std::endl;
    std::cout <<"result: " << result  << "\n" << std::endl;

    OrbitalTensor currentTensorCM((*itJPCLS)->L);
    currentTensorCM.SetP4(pbarCM4Vec,pCM4Vec);

    PolVector currentPolVecCM((*itJPCLS)->J);
    currentPolVecCM.SetP4(allCM4Vec, allCM4Vec.M());

    Tensor<complex<double> > resultCM= currentPolVecCM(0).Conjugate() | currentTensorCM;
    resultCM*= Psi0_m12_12CM;
    std::cout <<"resultCM: " << resultCM  << "\n" << std::endl;
  }

  std::cout <<"\n\n\n" << "*** triplet 0 states ******" <<"\n\n\n" << std::endl;
  std::vector< std::shared_ptr<const JPCLS> > jpclsStatesTriplet0=thepbarpState.triplet0_JPCLS_States();
  for(itJPCLS=jpclsStatesTriplet0.begin(); itJPCLS!=jpclsStatesTriplet0.end(); ++itJPCLS){
    OrbitalTensor currentTensor((*itJPCLS)->L);
    currentTensor.SetP4(pbar4Vec,p4Vec);

    PolVector currentPolVec((*itJPCLS)->J);
    currentPolVec.SetP4(allVec, allVec.M());

    OrbitalTensor currentTensorCM((*itJPCLS)->L);
    currentTensorCM.SetP4(pbarCM4Vec,pCM4Vec);

    PolVector currentPolVecCM((*itJPCLS)->J);
    currentPolVecCM.SetP4(allCM4Vec, allCM4Vec.M());

    Tensor<complex<double> > result;
    if ((*itJPCLS)->J > (*itJPCLS)->L) result= currentPolVec(0).Conjugate() | (Psi1_m12_12 % currentTensor);
    else result= currentPolVec(0).Conjugate() | (Psi1_m12_12 *  currentTensor );

    std::cout <<"\n" << (*itJPCLS)->name() << std::endl;
    std::cout <<"result Psi1 -1/2 1/2: " << result  << "\n" << std::endl;
 
    Tensor<complex<double> > result1;
    if ((*itJPCLS)->J > (*itJPCLS)->L) result1=currentPolVec(0).Conjugate() | (Psi1_12_m12 % currentTensor);
    else result1=currentPolVec(0).Conjugate() | (Psi1_12_m12 *  currentTensor );
    std::cout <<"result Psi1 1/2 -1/2: " << result1  << "\n" << std::endl;

    Tensor<complex<double> > resultCM;
    if ((*itJPCLS)->J > (*itJPCLS)->L) resultCM= currentPolVecCM(0).Conjugate() | (Psi1_m12_12CM % currentTensorCM);
    else resultCM= currentPolVecCM(0).Conjugate() | (Psi1_m12_12CM *  currentTensorCM );

    std::cout <<"resultCM Psi1 -1/2 1/2: " << resultCM  << "\n" << std::endl;
 
    Tensor<complex<double> > result1CM;
    if ((*itJPCLS)->J > (*itJPCLS)->L) result1CM= currentPolVecCM(0).Conjugate() | (Psi1_12_m12CM % currentTensorCM);
   else result1CM= currentPolVecCM(0).Conjugate() | (Psi1_12_m12CM *  currentTensorCM );
    std::cout <<"resultCM Psi1 1/2 -1/2: " << result1CM  << "\n" << std::endl; 
  }

  LeviCivitaTensor lctTensor;
  Tensor<complex<double> > levipTensor=lctTensor*allVec;

  Tensor<complex<double> > levipTensorCM=lctTensor*allCM4Vec;
  
 std::cout <<"\n\n\n" << "*** triplet 1 states ******" <<"\n\n\n" << std::endl;
  std::vector< std::shared_ptr<const JPCLS> > jpclsStatesTriplet1=thepbarpState.tripletp1_JPCLS_States();
  for(itJPCLS=jpclsStatesTriplet1.begin(); itJPCLS!=jpclsStatesTriplet1.end(); ++itJPCLS){
    OrbitalTensor currentTensor((*itJPCLS)->L);
    currentTensor.SetP4(pbar4Vec,p4Vec);

    PolVector currentPolVec((*itJPCLS)->J);
    currentPolVec.SetP4(allVec, allVec.M());

    OrbitalTensor currentTensorCM((*itJPCLS)->L);
    currentTensorCM.SetP4(pbarCM4Vec,pCM4Vec);

    PolVector currentPolVecCM((*itJPCLS)->J);
    currentPolVecCM.SetP4(allCM4Vec, allCM4Vec.M());

    Tensor<complex<double> > result;
    Tensor<complex<double> > resultm1m1;
    Tensor<complex<double> > resultCM;
    Tensor<complex<double> > resultm1m1CM;
    
    if ((*itJPCLS)->J > (*itJPCLS)->L){
      result= currentPolVec(1).Conjugate() | (Psi1_12_12 % currentTensor);
      resultm1m1= currentPolVec(-1).Conjugate() | (Psi1_m12_m12 % currentTensor);
      resultCM= currentPolVecCM(1).Conjugate() | (Psi1_12_12CM % currentTensorCM);
      resultm1m1CM= currentPolVecCM(-1).Conjugate() | (Psi1_m12_m12CM % currentTensorCM);
    }
    else if ((*itJPCLS)->J == (*itJPCLS)->L){
      Tensor<complex<double> > tmpResult=levipTensor*Psi1_12_12;
      Tensor<complex<double> > tmpResult1 = tmpResult.Permute(0, tmpResult.Rank()-1)*currentTensor;
      result= currentPolVec(1).Conjugate() | tmpResult1;

      Tensor<complex<double> > tmpResultm1m1=levipTensor*Psi1_m12_m12;
      Tensor<complex<double> > tmpResult1m1m1=tmpResultm1m1.Permute(0, tmpResultm1m1.Rank()-1)*currentTensor; 
      resultm1m1= currentPolVec(-1).Conjugate() | tmpResult1m1m1;

      Tensor<complex<double> > tmpResultCM=levipTensorCM*Psi1_12_12CM;
      Tensor<complex<double> > tmpResult1CM=tmpResultCM.Permute(0, tmpResultCM.Rank()-1)*currentTensorCM;
      resultCM= currentPolVecCM(1).Conjugate() | tmpResult1CM;

      Tensor<complex<double> > tmpResultm1m1CM=levipTensorCM*Psi1_m12_m12CM;
      Tensor<complex<double> > tmpResult1m1m1CM=tmpResultm1m1CM.Permute(0, tmpResultm1m1CM.Rank()-1)*currentTensorCM; 
      resultm1m1CM= currentPolVecCM(-1).Conjugate() | tmpResult1m1m1CM;      
    }
    else{
      result= currentPolVec(1).Conjugate() | (Psi1_12_12 *  currentTensor );
      resultm1m1= currentPolVec(-1).Conjugate() | (Psi1_m12_m12 *  currentTensor );
      resultCM= currentPolVecCM(1).Conjugate() | (Psi1_12_12CM *  currentTensorCM );
      resultm1m1CM= currentPolVecCM(-1).Conjugate() | (Psi1_m12_m12CM *  currentTensorCM );
    }

    std::cout <<"\n" << (*itJPCLS)->name() << std::endl;
    std::cout <<"result M1 Psi 1/2 1/2: " << result  << "\n" << std::endl;
    std::cout <<"result M-1 Psi -1/2 -1/2: " << resultm1m1  << "\n" << std::endl;
    std::cout <<"resultCM M1 Psi 1/2 1/2: " << resultCM  << "\n" << std::endl;
    std::cout <<"resultCM M-1 Psi -1/2 -1/2: " << resultm1m1CM  << "\n" << std::endl;
  }

  return 0;
}

