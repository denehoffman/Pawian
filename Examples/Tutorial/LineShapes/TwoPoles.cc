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
#include "Examples/Tutorial/LineShapes/TwoPoles.hh"
//#include "qft++/matrix/KpoleMatrix.hh"
//#include "qft++/matrix/TMatrix.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/KMatrixNonRel.hh"
#include "PwaDynamics/TMatrixNonRel.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceIsobar.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

TwoPoles::TwoPoles(double MassRes1, double Width1, double MassRes2, double Width2, double deltaMass) :
  _theTFile(0)
{

  if (MassRes1 < MassRes2){
    _massLow = MassRes1;
    _widthLow = Width1;
    _massHigh = MassRes2;
    _widthHigh = Width2;
    
  }
  else{
    _massLow = MassRes2;
    _widthLow = Width2;
    _massHigh = MassRes1;
    _widthHigh = Width1;
    
  }

  _theTFile=new TFile("TwoPoles.root","recreate");
 
  _massShapeHist= new TH1F("_massShapeHist", "BreitWigner",301, _massLow-deltaMass, _massHigh+deltaMass);
  _massShapeKmatrHist= new TH1F("_massShapeKmatrHist", "mass shape Kmatrix",301, _massLow-deltaMass, _massHigh+deltaMass);
  _massShapeKmatrHistRel= new TH1F("_massShapeKmatrHistRel", "mass shape rel. Kmatrix",301, _massLow-deltaMass, _massHigh+deltaMass);

  _argandHist=new TH2F("_argandHist","Argand plot",301, -1., 1., 301, 0., 1.3);
  _argandHist->SetXTitle("Re(Bw)");
  _argandHist->SetYTitle("Im(Bw)");
  _argandHist->SetMarkerStyle(6);
  
  _argandKmatrHist=new TH2F("_argandKmatrHist","Argand plot K matrix",301, -1., 1., 301, 0., 1.3);
  _argandKmatrHist->SetXTitle("Re(Bw)");
  _argandKmatrHist->SetYTitle("Im(Bw)");
  _argandKmatrHist->SetMarkerStyle(6);

  _argandKmatrHistRel=new TH2F("_argandKmatrHistRel","Argand plot rel. K matrix",301, -1., 1., 301, 0., 1.3);
  _argandKmatrHistRel->SetXTitle("Re(Bw)");
  _argandKmatrHistRel->SetYTitle("Im(Bw)");
  _argandKmatrHistRel->SetMarkerStyle(6);

  
  //initialize here the phase shift histogramms _phaseHistMap
  
  _phaseHist=new TH2F("_phaseHist", "phase shift",301, _massLow-deltaMass, _massHigh+deltaMass, 301, 0., 3.1415);
  _phaseHist->SetXTitle("mass");
  _phaseHist->SetYTitle("phase #delta(#phi)");
  _phaseHist->SetMarkerStyle(6);
  
  _phaseKmatrHist=new TH2F("_phaseKmatrHist", "phase shift Kmatrix",301, _massLow-deltaMass, _massHigh+deltaMass, 301, 0., 3.1415);
  _phaseKmatrHist->SetXTitle("mass");
  _phaseKmatrHist->SetYTitle("phase #delta(#phi)");
  _phaseKmatrHist->SetMarkerStyle(6);
  
  _phaseKmatrHistRel=new TH2F("_phaseKmatrHistRel", "phase shift rel Kmatrix",301, _massLow-deltaMass, _massHigh+deltaMass, 301, 0., 3.1415);
  _phaseKmatrHistRel->SetXTitle("mass");
  _phaseKmatrHistRel->SetYTitle("phase #delta(#phi)");
  _phaseKmatrHistRel->SetMarkerStyle(6);
  // break up momentum q0

  double piMass=0.135;
//   double rhoLow=sqrt( ( _massLow*_massLow-(piMass+piMass)*(piMass+piMass) ) 
// 		    * ( _massLow*_massLow - (piMass-piMass) * (piMass-piMass) ) ) / (_massLow);
  
//   double rhoHigh=sqrt( ( _massHigh*_massHigh-(piMass+piMass)*(piMass+piMass) ) 
// 		     * ( _massHigh*_massHigh - (piMass-piMass) * (piMass-piMass) ) ) / (_massHigh);

  double rhoLow=sqrt( (1-(2.*piMass/_massLow)*(2.*piMass/_massLow)) ); 
  
  double rhoHigh=sqrt( (1-(2.*piMass/_massHigh)*(2.*piMass/_massHigh)) ); 

  double gLow=sqrt(_massLow*_widthLow/rhoLow);
  double gHigh=sqrt(_massHigh*_widthHigh/rhoHigh);  

  pair <double, double> theMassDecPair(piMass, piMass);
  vector <pair <double, double> > theMassDecPairVec;
  theMassDecPairVec.push_back(theMassDecPair);

  vector<double> low_gis;
  low_gis.push_back(gLow);
  // KpoleMatrix kPoleLow(low_gis, _massLow, theMassDecPairVec);
  std::shared_ptr<KPole> kPoleLow(new KPole(low_gis, _massLow));

  vector<double> high_gis;
  high_gis.push_back(gHigh);
  // KpoleMatrix kPoleHigh(high_gis, _massHigh, theMassDecPairVec);
  std::shared_ptr<KPole> kPoleHigh(new KPole(high_gis, _massHigh));

  // vector<KpoleMatrix> kPoles;
  // kPoles.push_back(kPoleLow);
  // kPoles.push_back(kPoleHigh);
  std::vector<std::shared_ptr<KPole> > kPoles;
   kPoles.push_back(kPoleLow);
   kPoles.push_back(kPoleHigh);

   std::vector<std::shared_ptr<AbsPhaseSpace> > phpVec;
   std::shared_ptr<AbsPhaseSpace> firstPhp(new PhaseSpaceIsobar(piMass, piMass));
   phpVec.push_back(firstPhp); 

   // TMatrix theTMatrix(kPoles);
   std::shared_ptr<KMatrixRel> theKMatrixRel( new KMatrixRel(kPoles, phpVec));
   TMatrixRel theTMatrixRel(theKMatrixRel);


   std::shared_ptr<KMatrixNonRel> theKMatrixNonRel( new KMatrixNonRel(kPoles, phpVec));
   TMatrixNonRel theTMatrixNonRel(theKMatrixNonRel);  

  double stepSize=( (_massHigh+deltaMass)-(_massLow-deltaMass) )/300.;  




  for (double massIt=_massLow-deltaMass; massIt<_massHigh+deltaMass; massIt+=stepSize){
    
    Vector4<double>  res4V(massIt, 0., 0., 0.);
    
    complex<double> currentBWs=BreitWigner(res4V, _massLow, _widthLow)+BreitWigner(res4V, _massHigh, _widthHigh);
    double weight=norm(currentBWs);
    _massShapeHist->Fill(massIt,weight);
    
    _argandHist->Fill(currentBWs.real(),currentBWs.imag());
    
    //fill here the phase shift histogramms: hint "atan2(imag,real)"
    
    _phaseHist->Fill(massIt, atan2(currentBWs.imag(), currentBWs.real()));

    //Kmatrix
  // break up momentum q0
 //    double piMass=0.135;
//     double qlow=sqrt( ( _massLow*_massLow-(piMass+piMass)*(piMass+piMass) ) 
// 		      * ( _massLow*_massLow - (piMass-piMass) * (piMass-piMass) ) ) / (_massLow);

//     double qhigh=sqrt( ( _massHigh*_massHigh-(piMass+piMass)*(piMass+piMass) ) 
// 		       * ( _massHigh*_massHigh - (piMass-piMass) * (piMass-piMass) ) ) / (_massHigh);

    // break up momentum qab
 
    // phase space factor rhoLow
//     double  rhoLow=2*qlow/_massLow;

//    double  rhoHigh=2*qhigh/_massHigh;
    // phase space factor rhoMab
//     double  rhoMab=2*qab/massIt;


    theTMatrixNonRel.evalMatrix(massIt);
    complex<double> currentTmatrix=theTMatrixNonRel(0,0);
    cout << "TMatrix:\t" <<  currentTmatrix << endl;
    weight=norm(currentTmatrix); 
    _massShapeKmatrHist->Fill(massIt,weight);
    _argandKmatrHist->Fill(currentTmatrix.real(),currentTmatrix.imag()); 
    _phaseKmatrHist->Fill(massIt, atan2(currentTmatrix.imag(), currentTmatrix.real()));

    theTMatrixRel.evalMatrix(massIt);
    complex<double> currentTmatrixRel=theTMatrixRel(0,0);
    cout << "TMatrixRel:\t" <<  currentTmatrixRel << endl;
    weight=norm(currentTmatrixRel); 
    _massShapeKmatrHistRel->Fill(massIt,weight);
    _argandKmatrHistRel->Fill(currentTmatrixRel.real(),currentTmatrixRel.imag()); 
    _phaseKmatrHistRel->Fill(massIt, atan2(currentTmatrixRel.imag(), currentTmatrixRel.real()));
    
  }


}

TwoPoles::~TwoPoles()
{
   _theTFile->Write();
   _theTFile->Close();
}

