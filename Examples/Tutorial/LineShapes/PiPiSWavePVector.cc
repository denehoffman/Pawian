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
#include <boost/multi_array.hpp>
#include "Examples/Tutorial/LineShapes/PiPiSWavePVector.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/TMatrixBase.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/TMatrixNonRel.hh"
#include "PwaDynamics/KPole.hh"
//#include "PwaDynamics/KPoleNonRel.hh"
#include "PwaDynamics/KMatrixSlowAdlerCorRel.hh"
#include "PwaDynamics/KMatrixNonRel.hh"

#include "PwaDynamics/PPole.hh"
#include "PwaDynamics/PVectorRel.hh"
#include "PwaDynamics/PVectorSlowCorRel.hh"
#include "PwaDynamics/FVector.hh"
#include "PwaDynamics/FVectorPiPiS.hh"


#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

PiPiSWavePVector::PiPiSWavePVector() :
  _theTFile(0)
{
  std::string rootFileName="./PiPiSWavePVector.root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");

  std::vector< complex<double> > fProd;
  fProd.resize(5);
  // fProd[0]=complex <double>(- 2.19449   , -7.62666);
  // fProd[1]=complex <double>(1.87981, -.628378)*fProd[0];
  // fProd[2]=complex <double>(4.3242 , 2.75019)*fProd[0];
  // fProd[3]=complex <double>(3.22336,  .271048)*fProd[0];
  // fProd[4]=complex <double>(0.,0.);

  // fProd[0]=complex <double>(-0.798456,-0.327471);
  // fProd[1]=complex <double>(0.,0.);
  // fProd[2]=complex <double>(-2.8991,2.40688);
  // fProd[3]=complex <double>(-0.622515,2.08818);
  // fProd[4]=complex <double>(-1.35335,-2.04469);

  // fProd[0]=complex <double>(-2.33219,1.06776);
  // fProd[1]=complex <double>(-5.93134,-2.15883);
  // fProd[2]=complex <double>(0.,0.);
  // fProd[3]=complex <double>(0.,0.);
  // fProd[4]=complex <double>(0.,0.);

  fProd[0]=complex <double>(0.,0.);
  fProd[1]=complex <double>(0.,0.);
  fProd[2]=complex <double>(0.,0.);
  fProd[3]=complex <double>(0.,0.);
  fProd[4]=complex <double>(0.,0.);

  std::vector<complex <double> > betaPars;
  betaPars.push_back(complex <double>(-0.798456,-0.327471));
  betaPars.push_back(complex <double>(-10.1735,-3.88488));
  betaPars.push_back(complex <double>(-23.6712,5.03146));
  betaPars.push_back(complex <double>(-0.0799351,9.15965));
  betaPars.push_back(complex <double>(0.,0.));

  // betaPars.push_back(complex <double>(1.,0.));
  // betaPars.push_back(complex <double>(-0.86784,-0.442187));
  // betaPars.push_back(complex <double>(0.0566089,0.752875));
  // betaPars.push_back(complex <double>(-0.85595,1.53784));
  // betaPars.push_back(complex <double>(0.,0.));

  // betaPars.push_back(complex <double>(1.,0.));
  // betaPars.push_back(complex <double>(0.32253,2.44986));
  // betaPars.push_back(complex <double>(-0.240159,1.05706));
  // betaPars.push_back(complex <double>(0.,0.));
  // betaPars.push_back(complex <double>(0.,0.));

  // betaPars.push_back(complex <double>(1.8223,-9.11972));
  // betaPars.push_back(complex <double>(-10.1735,-3.88488));
  // betaPars.push_back(complex <double>(-23.6712,-5.03146));
  // betaPars.push_back(complex <double>(-0.0799351,9.15965));
  // betaPars.push_back(complex <double>(0.,0.));


  int size=2000;
  double massMin=2.*0.135;
  double massMax=1.9;

  // double massMax=1.4;

  double stepSize=(massMax-massMin)/size;

  _absF0H1= new TH1F("_absF0H1","|#rho_{00}#hat{F}_{0}|",size+1, massMin, massMax);
  _absF0H1->SetYTitle("|#rho_{00}#hat{F}_{0}|");
  _sqrF0H1= new TH1F("_sqrF0H1","|#rho_{00}#hat{F}_{0}|^{2}",size+1, massMin, massMax);
  _sqrF0H1->SetYTitle("|#rho_{00}#hat{F}_{0}|^{2}");
  _argandH2=new TH2F("_argandH2","Argand plot",200, -10., 10., 200, -10., 10.);
  _argandH2->SetYTitle("Im(F_{0})");
  _argandH2->SetXTitle("Re(F_{0})");
  _phaseShiftH2=new TH2F("_phaseShiftH2", "phase shift",301, massMin, massMax, 301, 0., 3.1415);


  // //  double s0Prod=-3.;
  // std::shared_ptr<PVectorRel> thePVector(new PVectorSlowCorRel(pPoles, phpVecs, fProdMatr, s0Prod));

  // //std::shared_ptr<PVectorRel> thePVector(new PVectorRel(pPoles, phpVecs));

  std::shared_ptr<FVectorPiPiS> theFVector(new FVectorPiPiS());

  std::vector<complex <double> >::const_iterator itBeta;
  for( int i=0; i< int(betaPars.size()); ++i){
    theFVector->updateBeta(i, betaPars[i]);    
  } 

  std::vector<complex <double> >::const_iterator itFprod;
  for( int i=0; i< int(fProd.size()); ++i){
    theFVector->updateFprod(i, fProd[i]);    
  } 

  double s0Prod=-0.0737;
  theFVector->updateS0prod(s0Prod);

  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=theFVector->kMatrix()->phaseSpaceVec();  
  
  for (double mass=massMin; mass<massMax; mass+=stepSize){
    Vector4<double> mass4Vec(mass, 0.,0.,0.);
    
    theFVector->evalMatrix(mass);
    complex<double> currentValFvec=(*theFVector)(0,0);
    _absF0H1->Fill(mass4Vec.M(), sqrt(norm(currentValFvec)));
    _sqrF0H1->Fill(mass4Vec.M(), norm(currentValFvec) );
    _argandH2->Fill( currentValFvec.real(), currentValFvec.imag() );
    _phaseShiftH2->Fill(mass, atan2(currentValFvec.imag(), currentValFvec.real()));
  }
}

PiPiSWavePVector::~PiPiSWavePVector()
{
  _theTFile->Write();
  _theTFile->Close();
}

