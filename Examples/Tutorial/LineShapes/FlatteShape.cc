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
#include "Examples/Tutorial/LineShapes/FlatteShape.hh"
#include "PwaDynamics/Flatte.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceIsobar.hh"
//#include "qft++/matrix/KpoleMatrix.hh"
//#include "qft++/matrix/TMatrix.hh"
#include "qft++Extension/PawianUtils.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

FlatteShape::FlatteShape(std::string ptype, double g1, double g2) :
  _theTFile(0)
  ,_histShapeLow(0)
  ,_histShapeHigh(0)
  ,_pType(ptype)
{
  if(_pType!="a0" && _pType!="f0"){
    Alert <<"only a0 and f0 are supported" << endmsg;
    exit(1);
  }
  if(g1<0.  || g2<0.){
    Alert <<"requirement: g1>0. and g2>0." << endmsg;
    exit(1);
  }

  std::string rootFileName="./FlatteShape"+_pType+".root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");


  double mass1=0.;
  double mass2=0.;

  if(_pType=="a0"){
    mass1 = 0.1349766;
    mass2 = 0.547853;
  }
  else if(_pType=="f0"){
    mass1 = 0.1349766;
    mass2= 0.1349766;
  }

  std::pair <const double, const double> decPairLow=make_pair(mass1, mass2); 
  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;
  std::pair <const double, const double> decPairHigh=make_pair(KplusMass,K0Mass);


  int size=2000;
  double massMin=.5;
  double massMax=1.85;

//   double massMin=1.;
//   double massMax=1.004;
  double stepSize=(massMax-massMin)/size;

  _histShapeLow= new TH1F("_histShapeLow","hist low",size+1, massMin, massMax);
  _histShapeHigh= new TH1F("_histShapeHigh","hist high",size+1, massMin, massMax);
  _histShapeLowFlatteNew= new TH1F("_histShapeLowFlatteNew","hist high flatte new",size+1, massMin, massMax);
  _histShapeHighFlatteNew= new TH1F("_histShapeHighFlatteNew","hist high",size+1, massMin, massMax);
  _histShapeLowTest= new TH1F("_histShapeLowTest","hist low test",size+1, massMin, massMax);
  _histShapeHighTest= new TH1F("_histShapeHighTest","hist high test",size+1, massMin, massMax);
  _histShapeLowKmatr= new TH1F("_histShapeLowKmatr","hist low K-matrix",size+1, massMin, massMax);
  _histShapeHighKmatr= new TH1F("_histShapeHighKmatr","hist high K-matrix",size+1, massMin, massMax);
  _histShapeLowKmatrRel= new TH1F("_histShapeLowKmatrRel","hist low K-matrix relativistic",size+1, massMin, massMax);
  _histShapeHighKmatrRel= new TH1F("_histShapeHighKmatrRel","hist high K-matrix relativistic",size+1, massMin, massMax);

  _argandKmatrFlatteLowHist=new TH2F("_argandKmatrFlatteLowHist","Argand plot K matrix",301, -1., 1., 301, 0., 1.3);
  _argandKmatrFlatteLowHist->SetXTitle("Re(Bw)");
  _argandKmatrFlatteLowHist->SetYTitle("Im(Bw)");
  _argandKmatrFlatteLowHist->SetMarkerStyle(6);

  _argandKmatrLowHist=new TH2F("_argandKmatrLowHist","Argand plot K matrix",301, -1., 1., 301, 0., 1.3);
  _argandKmatrLowHist->SetXTitle("Re(Bw)");
  _argandKmatrLowHist->SetYTitle("Im(Bw)");
  _argandKmatrLowHist->SetMarkerStyle(6);

  _argandKmatrHighHist=new TH2F("_argandKmatrHighHist","Argand plot K matrix high",301, -1., 1., 301, 0., 1.3);
  _argandKmatrHighHist->SetXTitle("Re(Bw)");
  _argandKmatrHighHist->SetYTitle("Im(Bw)");
  _argandKmatrHighHist->SetMarkerStyle(6);

  _argandKmatrLowRelHist=new TH2F("_argandKmatrLowRelHist","Argand plot K matrix relativistic high",301, -1., 1., 301, 0., 1.3);
  _argandKmatrLowRelHist->SetXTitle("Re(Bw)");
  _argandKmatrLowRelHist->SetYTitle("Im(Bw)");
  _argandKmatrLowRelHist->SetMarkerStyle(6);

  _argandKmatrHighRelHist=new TH2F("_argandKmatrHighRelHist","Argand plot K matrix relativistic",301, -1., 1., 301, 0., 1.3);
  _argandKmatrHighRelHist->SetXTitle("Re(Bw)");
  _argandKmatrHighRelHist->SetYTitle("Im(Bw)");
  _argandKmatrHighRelHist->SetMarkerStyle(6);

  vector<double> the_gis;

  the_gis.push_back(g1);
  the_gis.push_back(g2);

  vector <pair <double, double> > theMassDecPairVec;
  theMassDecPairVec.push_back(decPairLow);
  theMassDecPairVec.push_back(decPairHigh);

  //    KpoleMatrix kPole(the_gis, 0.98, theMassDecPairVec);
  std::shared_ptr<KPole> kPole(new KPole(the_gis, 0.98));

   // vector<KpoleMatrix> kPoles;
   // kPoles.push_back(kPole);

  std::vector<std::shared_ptr<KPole> > kPoles;
   kPoles.push_back(kPole);

   std::vector<std::shared_ptr<AbsPhaseSpace> > phpVec;
   std::shared_ptr<AbsPhaseSpace> lowPhp(new PhaseSpaceIsobar(decPairLow.first, decPairLow.second));
   phpVec.push_back(lowPhp);
   std::shared_ptr<AbsPhaseSpace> kkPhp(new PhaseSpaceIsobar(KplusMass, K0Mass));
   phpVec.push_back(kkPhp);

   std::shared_ptr<KMatrixRel> theKMatrixRel( new KMatrixRel(kPoles, phpVec));


   TMatrixRel theTMatrix(theKMatrixRel);

   Flatte theFlatte(decPairLow, decPairHigh); 
  
   for (double mass=massMin; mass<massMax; mass+=stepSize){
     Vector4<double> mass4Vec(mass, 0.,0.,0.);
     complex<double> flatteLow=PawianQFT::FlatteFkt(mass4Vec, decPairLow, decPairHigh, 0.98, g1, g2);
     cout << "flatteLow" << flatteLow << "  norm: " << norm(flatteLow) << endl;
     _histShapeLow->Fill(mass4Vec.M(), norm(flatteLow));
     
     _argandKmatrFlatteLowHist->Fill(flatteLow.real(),flatteLow.imag());
    
     complex<double> flatteHigh=PawianQFT::FlatteFkt(mass4Vec, decPairHigh, decPairLow, 0.98, g2, g1);

     cout << "flatteHigh" << flatteHigh << "  norm: " << norm(flatteHigh) << endl;
     _histShapeHigh->Fill(mass4Vec.M(), norm(flatteHigh));

     complex<double> flatteLowNew=theFlatte.calcFirstChannel(mass, 0.98, g1, g2);
     _histShapeLowFlatteNew->Fill(mass4Vec.M(), norm(flatteLowNew));

     complex<double> flatteHighNew=theFlatte.calcSecondChannel(mass, 0.98, g1, g2);
     _histShapeHighFlatteNew->Fill(mass4Vec.M(), norm(flatteHighNew));
     
     theTMatrix.evalMatrix(mass);
     cout << "current TMatrix:\n " <<  theTMatrix << endl;
     
     complex<double> currentValLow=theTMatrix(0,0);
     cout << "currentMass= " << mass4Vec.M() << endl;
     cout << "currentValLow= " << currentValLow << "  norm: " <<  norm(currentValLow) << endl;
     complex<double> currentValHigh=theTMatrix(0,1); 
     cout << "currentValHigh= " << currentValHigh << "  norm: " <<  norm(currentValHigh) << endl;
 
     _histShapeLowKmatr->Fill(mass4Vec.M(), norm(currentValLow) );
     _argandKmatrLowHist->Fill(currentValLow.real(),currentValLow.imag());
     _histShapeHighKmatr->Fill(mass4Vec.M(), norm(currentValHigh) );
     _argandKmatrHighHist->Fill(currentValHigh.real(),currentValHigh.imag());

     theTMatrix.evalMatrix(mass);
     cout << "current TMatrixRel:\n " <<  theTMatrix << endl;

     complex<double> currentValLowRel=theTMatrix(0,0);
     cout << "currentValLowRel= " << currentValLowRel << "  norm: " <<  norm(currentValLowRel) << endl;
     complex<double> currentValHighRel=theTMatrix(0,1); 
     cout << "currentValHighRel= " << currentValHighRel << "  norm: " <<  norm(currentValHighRel) << endl;

     //     vector< complex<double> > rhoFactors=theTMatrix.currentRhoFactors(); 
     // complex<double> rhoDiv=rhoFactors[1]/rhoFactors[0];    
     vector< complex<double> > rhoFactors;
     rhoFactors.push_back(lowPhp->factor(mass4Vec.M()));
     rhoFactors.push_back(kkPhp->factor(mass4Vec.M()));
     
     _histShapeLowKmatrRel->Fill(mass4Vec.M(), norm(currentValLowRel)*norm(rhoFactors[0]) );
     _argandKmatrLowRelHist->Fill(sqrt(norm(rhoFactors[0]))*currentValLowRel.real(),sqrt(norm(rhoFactors[0]))*currentValLowRel.imag());
     _histShapeHighKmatrRel->Fill(mass4Vec.M(), norm(currentValHighRel)*norm(sqrt(rhoFactors[0]*rhoFactors[1])) );
     _argandKmatrHighRelHist->Fill(currentValHighRel.real()*sqrt(abs(rhoFactors[0]*rhoFactors[1])),currentValHighRel.imag()*sqrt(abs(rhoFactors[0]*rhoFactors[1])));

     complex<double> rhoPoleFactor1=PawianQFT::phaseSpaceFacDefault(mass4Vec.M(),decPairLow.first, decPairLow.second); 
     complex<double> rhoPoleFactor2=PawianQFT::phaseSpaceFacDefault(mass4Vec.M(),decPairHigh.first, decPairHigh.second);

     complex<double> imag(0.,1.);
     complex<double> denom=0.98*0.98-mass4Vec.M()*mass4Vec.M()-imag*(rhoPoleFactor1*g1*g1+rhoPoleFactor2*g2*g2);
     complex<double> T00= g1*g1/denom;
     complex<double> T01= g1*g2/denom;
     complex<double> T11= g2*g2/denom;   
     cout << "T00= " << T00 << "\t norm= " << norm(T00) << endl;
     cout << "T11= " << T11 << "\t norm= " << norm(T11) << endl;
     _histShapeLowTest->Fill(mass4Vec.M(), norm(T00) );
     _histShapeHighTest->Fill(mass4Vec.M(), norm(T01) ); 
     

//     complex<double> nom1=0.98*0.98-mass4Vec.M()*mass4Vec.M()+imag*(rhoPoleFactor1*g1*g1+rhoPoleFactor2*g2*g2);
//     double denom1=(0.98*0.98-mass4Vec.M()*mass4Vec.M())*(0.98*0.98-mass4Vec.M()*mass4Vec.M())+norm(rhoPoleFactor1*g1*g1+rhoPoleFactor2*g2*g2);
//      complex<double> T00= g1*g1*nom1/denom1;
//      complex<double> T01= g1*g2*nom1/denom1;
//      complex<double> T11= g2*g2*nom1/denom1;
//      cout << "T00= " << T00 << "\t norm= " << norm(T00) << endl;
//      cout << "T11= " << T11 << "\t norm= " << norm(T11) << endl;
//      _histShapeLowTest->Fill(mass4Vec.M(), norm(T00) );
//      _histShapeHighTest->Fill(mass4Vec.M(), norm(T01) );  

//      cout << "2*i*2= " << 2.*imag*2. << endl;
//      cout << "4*i  = " << 4.*imag << endl;     
   }
   
}

FlatteShape::~FlatteShape()
{
  _theTFile->Write();
  _theTFile->Close();
}

