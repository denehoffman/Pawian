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

#include "PwaDynamics/KMatrixRelBg.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"

KMatrixRelBg::KMatrixRelBg(vector<std::shared_ptr<KPole> > Kpoles, vector<std::shared_ptr<AbsPhaseSpace> > phpVecs, unsigned int orderBg, bool withAdler) :
  KMatrixRel(Kpoles, phpVecs)
  ,_withAdler(withAdler)
 {
   _orderBg=orderBg;
   _bgTerms.resize(_orderBg+1);
   for(unsigned int i=0; i<=_orderBg; ++i){
     _bgTerms.at(i).resize(phpVecs.size());
     for(unsigned int j=0; j<phpVecs.size(); ++j){
       _bgTerms.at(i).at(j).resize(phpVecs.size());
       for(unsigned int k=0; k<phpVecs.size(); ++k){
	 _bgTerms.at(i).at(j).at(k)=0.;
       }
     }
   }
 }


KMatrixRelBg::~KMatrixRelBg(){
}

void KMatrixRelBg::evalMatrix(const double mass){
   evalMatrixTemplate(mass);
}

void KMatrixRelBg::evalMatrix(const complex<double> mass){
   evalMatrixTemplate(mass);
}

template<typename MassType>
void KMatrixRelBg::evalMatrixTemplate(const MassType mass){

  MassType adlerTerm=1.;
  MassType s_hat=mass*mass;
  if(_withAdler){
    adlerTerm=(mass*mass-_s0Adler)/_snormAdler;
    s_hat=mass*mass/_snormAdler-1.;
  }

  Matrix< complex<double> > theKMatrix(NumRows(), NumRows());
  vector<std::shared_ptr<KPole> >::iterator it;
  for (it =_KPoles.begin(); it != _KPoles.end(); ++it){
    (*it)->evalMatrix(mass);
    theKMatrix += *(*it);
  }

  for (int i=0; i<NumRows(); ++i){
    for (int j=i; j<NumRows(); ++j){
      complex<double> currentBg(0.,0.);
      for (unsigned int k=0; k<=_orderBg; ++k){
	currentBg+=complex<double>(_bgTerms.at(k).at(i).at(j)*pow(s_hat,k));
      }
      this->operator()(i,j)=theKMatrix(i,j)+currentBg;
      this->operator()(i,j)*=adlerTerm;
      this->operator()(j,i)=this->operator()(i,j);
    }
  }
}

template void KMatrixRelBg::evalMatrixTemplate(const double mass);
template void KMatrixRelBg::evalMatrixTemplate(const complex<double> mass);
