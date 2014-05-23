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
#include "Examples/Tutorial/LineShapes/KPiSWaveTMatrix.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh" 
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/TMatrixBase.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/TMatrixNonRel.hh"
#include "PwaDynamics/KMatrixBase.hh"
#include "PwaDynamics/KMatrixKPiSFocus.hh"
#include "PwaDynamics/KPole.hh"
//#include "PwaDynamics/KPoleNonRel.hh"
#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/LASS.hh"


#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"


#include "ErrLogger/ErrLogger.hh"

KPiSWaveTMatrix::KPiSWaveTMatrix() :
  _theTFile(0)
{
  std::string rootFileName="./KPiSWaveTMatrix.root";
  _theTFile=new TFile(rootFileName.c_str(),"recreate");

  double mLASS=1.435;
  double gammaLASS=0.279;
  // double aLASS=1.9/1000.;
  // double rLASS=1.76/1000.;
  double aLASS=1.9;
  double rLASS=1.76;
  double BLASS=1.0;
  double phiB=0.;
  double RLASS=1.0;
  double phiR=0.;

  int size=2000;
  double massMin=0.135+0.5;
  double massMax=1.8;

  double stepSize=(massMax-massMin)/size;

  _KPiAmpRealH1= new TH1F("_KPiAmpRealH1","K #pi amp Re",size+1, massMin, massMax);
  _KPiAmpRealH1->SetYTitle("K #pi amp Re");
  _KPiAmpImagH1= new TH1F("_KPiAmpImagH1","K #pi amp Im",size+1, massMin, massMax);
  _KPiAmpImagH1->SetYTitle("K #pi amp Im");

  _KPiAmpRealLASSH1= new TH1F("_KPiAmpRealLASSH1","K #pi amp Re LASS",size+1, massMin, massMax);
  _KPiAmpRealLASSH1->SetYTitle("K #pi amp Re");
  _KPiAmpImagLASSH1= new TH1F("_KPiAmpImagLASSH1","K #pi amp Im LASS",size+1, massMin, massMax);
  _KPiAmpImagLASSH1->SetYTitle("K #pi amp Im");

  std::shared_ptr<KMatrixBase> theKMatrixIso12(new KMatrixKPiSFocus(1));
  std::shared_ptr<KMatrixBase> theKMatrixIso32(new KMatrixKPiSFocus(3));
  // std::shared_ptr<TMatrixRel> theTMatrix(new TMatrixRel(theKMatrix));
  std::shared_ptr<TMatrixRel> theTMatrix12(new TMatrixRel(theKMatrixIso12));
  std::shared_ptr<TMatrixRel> theTMatrix32(new TMatrixRel(theKMatrixIso32));

  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=theKMatrixIso12->phaseSpaceVec();

  
  for (double mass=massMin; mass<massMax; mass+=stepSize){
    Vector4<double> mass4Vec(mass, 0.,0.,0.);

    complex<double> currentLASS=LASS::K0star_1430(mass, mLASS, gammaLASS, aLASS, rLASS, BLASS, phiB, RLASS, phiR);
    _KPiAmpRealLASSH1->Fill(mass, currentLASS.real());
    _KPiAmpImagLASSH1->Fill(mass, currentLASS.imag());
     
    theTMatrix12->evalMatrix(mass);
    theTMatrix32->evalMatrix(mass);        

    complex<double> currentVal=(*theTMatrix12)(0,0)+0.5*(*theTMatrix32)(0,0);
    //    complex<double> currentVal=(*theTMatrix12)(0,0);
       _KPiAmpRealH1->Fill(mass, sqrt(thePhpVecs[0]->factor(mass4Vec.M()).real())*currentVal.real());    
       _KPiAmpImagH1->Fill(mass, sqrt(thePhpVecs[0]->factor(mass4Vec.M()).real())*currentVal.imag());
  }
}

KPiSWaveTMatrix::~KPiSWaveTMatrix()
{
  _theTFile->Write();
  _theTFile->Close();
}

