//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#include <boost/math/distributions/normal.hpp> // for normal_distribution
#include "PwaDynamics/KMatrixFunctions.hh"
#include "ErrLogger/ErrLogger.hh"
#include "qft++Extension/PawianUtils.hh"

#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/KMatrixRelBg.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "Utils/PawianConstants.hh"

double KMatrixFunctions::twoDeltaArgand(std::shared_ptr<TMatrixRel> tMatr, unsigned int projectionIndex, double mass){
 //  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=tMatr->kMatrix()->phaseSpaceVec();
 // tMatr->evalMatrix(mass); 
 // complex<double> currentTiiRel_rho=(*tMatr)(projectionIndex, projectionIndex)*thePhpVecs.at(projectionIndex)->factor(mass).real();
 // double currentReERel = currentTiiRel_rho.real();
 // double currentImERel = currentTiiRel_rho.imag() - 0.5;

 // double twodeltaRel = atan2(currentImERel, fabs(currentReERel))*PawianConstants::radToDeg + 90.0;
 // if (currentTiiRel_rho.real()  < 0.0) {twodeltaRel = 360.0 - twodeltaRel;}
 // // unsigned int noRots= calcNoRots(tMatr, projectionIndex, mass);
 // // twodeltaRel+=noRots*360.; 
 // return twodeltaRel;

return KMatrixFunctions::deltaArgand(tMatr, projectionIndex, mass)*2.; 
}

double KMatrixFunctions::deltaArgand(std::shared_ptr<TMatrixRel> tMatr, unsigned int projectionIndex, double mass){
  // return  KMatrixFunctions::twoDeltaArgand(tMatr, projectionIndex, mass)/2.;
  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=tMatr->kMatrix()->phaseSpaceVec();
  tMatr->evalMatrix(mass); 
  complex<double> currentTiiRel_rho=(*tMatr)(projectionIndex, projectionIndex)*thePhpVecs.at(projectionIndex)->factor(mass).real();
  double currentReERel = currentTiiRel_rho.real();
  double currentImERel = currentTiiRel_rho.imag() - 0.5;
  
  double deltaRel = 0.5*atan2(currentImERel, fabs(currentReERel))*PawianConstants::radToDeg + 45.0;
  if (currentTiiRel_rho.real()  < 0.0) {deltaRel = 180.0 - deltaRel;}
  return deltaRel;
}

double KMatrixFunctions::deltaArgandWSigma(std::shared_ptr<TMatrixRel> tMatr, unsigned int projectionIndex, double mass, double sigma, unsigned int noOfPoints){
  boost::math::normal theGaus(mass, sigma); //gaussian
  //calculate delta for a couple of points within 2 sigma around the mean value and 
  //avarage it according to an gaussian
  double currentPdf=boost::math::pdf(theGaus, mass);
  double deltaMean=deltaArgand(tMatr, projectionIndex, mass);
  if(noOfPoints==0) return deltaMean;

  double deltaAv=deltaMean*currentPdf;
  double normalisation=currentPdf;

  // InfoMsg << "\nmean Mass: " << mass << "\tcurrentPdf: " << currentPdf << "\tdeltaArgand: " << deltaMean << endmsg;

  double stepSize=2.*sigma/(2.*noOfPoints);

  for(double currentMass = mass-2.*sigma; currentMass < mass+2.*sigma+1.e-8; currentMass+=stepSize){
    double currentDelta=deltaArgand(tMatr, projectionIndex, currentMass);
    if( (deltaMean-currentDelta) > 90. ) currentDelta+=180.;
    else if( (deltaMean-currentDelta) < -90. ) currentDelta-=180; 
    currentPdf=boost::math::pdf(theGaus, currentMass);
    deltaAv+=deltaArgand(tMatr, projectionIndex, currentMass)*currentPdf;     
    normalisation+=currentPdf;
    // InfoMsg << "currentMass: " << currentMass << "\tcurrentPdf: " << currentPdf << "\tdeltaArgand: " << currentDelta << endmsg; 
    
  }

  deltaAv=deltaAv/normalisation;  
  return deltaAv;   
}

unsigned int KMatrixFunctions::calcNoRots(std::shared_ptr<TMatrixRel> tMatr, unsigned int projectionIndex, double mass){
  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=tMatr->kMatrix()->phaseSpaceVec();
  std::shared_ptr<AbsPhaseSpace> thePhp=thePhpVecs.at(projectionIndex);
  double threshold= thePhp->thresholdMass();

 int noOfLoops=0;
 double stepSize=0.03;
 double lastCurrentImERel=-10.;
 double currentImERel=-10.;
 for(double currentMass=threshold; currentMass < mass+1.1*stepSize; currentMass+=stepSize){
   tMatr->evalMatrix(currentMass);
   complex<double> currentTiiRel_rho=(*tMatr)(projectionIndex, projectionIndex)*thePhpVecs.at(projectionIndex)->factor(mass).real();
   double nextImERel = currentTiiRel_rho.imag();
   if(nextImERel>currentImERel && lastCurrentImERel>currentImERel) noOfLoops++;
   lastCurrentImERel=currentImERel;
   currentImERel=nextImERel;
 }
 return noOfLoops; 
}

unsigned int KMatrixFunctions::noOfPhaseRotationsArgand(std::shared_ptr<TMatrixRel> tMatr, unsigned int projectionIndex, double mass){
  vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=tMatr->kMatrix()->phaseSpaceVec();
  std::shared_ptr<AbsPhaseSpace> thePhp=thePhpVecs.at(projectionIndex);
   double threshold= thePhp->thresholdMass();

  //now loop from threshold to the mass of interest
  int noOfLoops=0;
  double oldReERel=1.;
  double stepSize=0.01;
  for(double currentMass=threshold; currentMass < mass+stepSize; currentMass+=stepSize){
    tMatr->evalMatrix(currentMass);
    complex<double> currentTiiRel_rho=(*tMatr)(projectionIndex, projectionIndex)*thePhp->factor(currentMass).real();
    double currentReERel = currentTiiRel_rho.real();
    if (oldReERel < 0.0 && currentReERel > 0.0) {noOfLoops++;}
    //    InfoMsg << "currentReERel: " << currentReERel << "\toldReERel: " << oldReERel << endmsg;
    oldReERel = currentReERel;   
  }
  // InfoMsg << "mass: " << mass << "\tnoOfLoops: " << noOfLoops << endmsg;  
  return noOfLoops; 
}

void KMatrixFunctions::fillRotationArgandMap(std::shared_ptr<TMatrixRel> tMatr, unsigned int projectionIndex, std::map<unsigned int, double>& toFill, double massMax){
 vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=tMatr->kMatrix()->phaseSpaceVec();
 std::shared_ptr<AbsPhaseSpace> thePhp=thePhpVecs.at(projectionIndex);
  double threshold= thePhp->thresholdMass();

 //now loop from threshold to the mass of interest
 int noOfLoops=0;
 double oldReERel=1.;
 double stepSize=0.001;
 for(double currentMass=threshold; currentMass < massMax+stepSize; currentMass+=stepSize){
   tMatr->evalMatrix(currentMass);
   complex<double> currentTiiRel_rho=(*tMatr)(projectionIndex, projectionIndex)*thePhp->factor(currentMass).real();
   double currentReERel = currentTiiRel_rho.real();
   if (oldReERel < 0.0 && currentReERel > 0.0) {
     noOfLoops++;
     toFill[noOfLoops]=currentMass;
   }
   //    InfoMsg << "currentReERel: " << currentReERel << "\toldReERel: " << oldReERel << endmsg;
   oldReERel = currentReERel;   
 }

 // InfoMsg << "the filled map contains:" << endmsg;
 // std::map<unsigned int, double >::iterator it;
 // for(it=toFill.begin(); it!=toFill.end(); ++it){
 //   InfoMsg << it->first << "\t" << it->second << endmsg;
 // }

}

// void KMatrixFunctions::fillRotationArgandMap(std::shared_ptr<TMatrixRel> tMatr, unsigned int projectionIndex, std::map<unsigned int, double>& toFill, double massMax){
//    vector<std::shared_ptr<AbsPhaseSpace> > thePhpVecs=tMatr->kMatrix()->phaseSpaceVec();
//    std::shared_ptr<AbsPhaseSpace> thePhp=thePhpVecs.at(projectionIndex);
//     double threshold= thePhp->thresholdMass();

//    //now loop from threshold to the mass of interest
//    int noOfLoops=0;
//    double oldPhase=0.;
//    double stepSize=0.001;
//    for(double currentMass=threshold; currentMass < massMax+stepSize; currentMass+=stepSize){
//      //     double currentPhase=KMatrixFunctions::twoDeltaArgand(tMatr, projectionIndex, currentMass);
//      double currentPhase=KMatrixFunctions::deltaArgand(tMatr, projectionIndex, currentMass);
//      //phase should be between 0. and 360 degree
//      while(currentPhase>180.) currentPhase-=180.;
//      while(currentPhase<0.) currentPhase+=180.;

//      //     if ( oldPhase>90. && (currentPhase-oldPhase) < -30.) {
//      if ( (currentPhase-oldPhase) < -10.) {
//        // InfoMsg << "currentPhase: " << currentPhase << "\toldPhase: " << oldPhase << endmsg;
//        noOfLoops++;
//        toFill[noOfLoops]=currentMass;
//      }
//      //    InfoMsg << "currentReERel: " << currentReERel << "\toldReERel: " << oldReERel << endmsg;
//      oldPhase = currentPhase;   
//    }

//   // InfoMsg << "the filled map contains:" << endmsg;
//   // std::map<unsigned int, double >::iterator it;
//   // for(it=toFill.begin(); it!=toFill.end(); ++it){
//   //   InfoMsg << it->first << "\t" << it->second << endmsg;
//   // } 
  
// }





