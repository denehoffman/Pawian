//************************************************************************//
//									  //
//  Copyright 2020 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#include "PwaDynamics/PVectorBgOmnes.hh"
#include "PwaDynamics/PPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"
#include "ErrLogger/ErrLogger.hh"

PVectorBgOmnes::PVectorBgOmnes(vector<std::shared_ptr<PPole> > Ppoles, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs, unsigned int orderBg) :
  PVectorRelBg(Ppoles, phpVecs, orderBg)
 {
 }

PVectorBgOmnes::~PVectorBgOmnes(){
}

void PVectorBgOmnes::evalMatrix(const double mass, Spin OrbMom){

  vector<std::shared_ptr<PPole> >::iterator it;
  for (it =_Ppoles.begin(); it != _Ppoles.end(); ++it){
    (*it)->evalMatrix(mass, OrbMom);
  }

  for (int i=0; i<NumRows(); ++i){
    complex<double> currentBg(0.,0.);
    for(unsigned int bgOrder=0; bgOrder<=_orderBgP; ++bgOrder) currentBg+=_bgPTerms.at(bgOrder).at(i)*pow(mass*mass,bgOrder);

    complex<double> currentP(0.,0.);
    for (it =_Ppoles.begin(); it != _Ppoles.end(); ++it){
      (*it)->evalMatrix(mass, OrbMom);
      currentP += (*(*it))(i,0);
    }
    currentP+currentBg;
    this->operator()(i,0)=currentP;
  }

}

void PVectorBgOmnes::evalMatrix(const complex<double> mass, Spin OrbMom){

  vector<std::shared_ptr<PPole> >::iterator it;
  for (it =_Ppoles.begin(); it != _Ppoles.end(); ++it){
    (*it)->evalMatrix(mass, OrbMom);
  }

  for (int i=0; i<NumRows(); ++i){
    complex<double> currentBg(0.,0.);
    for(unsigned int bgOrder=0; bgOrder<=_orderBgP; ++bgOrder) currentBg+=_bgPTerms.at(bgOrder).at(i)*pow(mass*mass,bgOrder);

    complex<double> currentP(0.,0.);
    for (it =_Ppoles.begin(); it != _Ppoles.end(); ++it){
      (*it)->evalMatrix(mass, OrbMom);
      currentP += ((*(*it))(i,0)+currentBg);
    }
    this->operator()(i,0)=currentP;
  }

}

