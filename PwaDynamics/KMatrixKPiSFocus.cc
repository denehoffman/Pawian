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

#include "PwaDynamics/KMatrixKPiSFocus.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceIsobar.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"
#include "ErrLogger/ErrLogger.hh"

KMatrixKPiSFocus::KMatrixKPiSFocus(int isospin2) :
  KMatrixBase(2, 2)
  ,_sAdler0(0.23)
  ,_sNorm(0.493677*0.49367+0.13957018*0.13957018)
{
  const double piMass=0.1349766;
  const double KplusMass=0.493677;
  const double etaprimeMass= 0.95778;

  std::shared_ptr<AbsPhaseSpace> kpiPhp(new PhaseSpaceIsobar(KplusMass, piMass));
  std::shared_ptr<AbsPhaseSpace> ketapPhp(new PhaseSpaceIsobar(KplusMass, etaprimeMass));
  _phpVecs.push_back(kpiPhp);
  _phpVecs.push_back(ketapPhp);

  _aScatPtr=std::shared_ptr< boost::multi_array<double, 2> >(new boost::multi_array<double, 2>(boost::extents[NumRows()][NumCols()]));
  _bScatPtr=std::shared_ptr< boost::multi_array<double, 2> >(new boost::multi_array<double, 2>(boost::extents[NumRows()][NumCols()]));
  _cScatPtr=std::shared_ptr< boost::multi_array<double, 2> >(new boost::multi_array<double, 2>(boost::extents[NumRows()][NumCols()]));

  for (int i=0; i<NumRows(); ++i){
    for (int j=0; j<NumCols(); ++j){
      (*_aScatPtr)[i][j]=0.;
      (*_bScatPtr)[i][j]=0.;
      (*_cScatPtr)[i][j]=0.;
    }
  }

  if (isospin2==1) init2IsoSpin1();
  else if (isospin2==3) init2IsoSpin3();
  else {
    Alert << "K pi S wave with isospin\t" << isospin2 << "/2 does not exist" << endmsg;
    exit(0);
  }
}

KMatrixKPiSFocus::~KMatrixKPiSFocus(){
}

void KMatrixKPiSFocus::evalMatrix(const double mass, Spin OrbMom){
   evalMatrixTemplate(mass);
}

void KMatrixKPiSFocus::evalMatrix(const complex<double> mass, Spin OrbMom){
   evalMatrixTemplate(mass);
}

template<typename MassType>
void KMatrixKPiSFocus::evalMatrixTemplate(const MassType mass, Spin OrbMom){
 
  Matrix< complex<double> > theKMatrix(NumRows(), NumRows());
  vector<std::shared_ptr<KPole> >::iterator it;
  for (it =_KPoles.begin(); it != _KPoles.end(); ++it){
    (*it)->evalMatrix(mass);
    theKMatrix += *(*it);
  }

  MassType s_hat=mass*mass/(_sNorm)-1.;  
  
  for (int i=0; i<theKMatrix.NumRows(); ++i){
    for (int j=0; j<theKMatrix.NumCols(); ++j){
      this->operator()(i,j)=theKMatrix(i,j)+complex<double> ( (*_aScatPtr)[i][j] + (*_bScatPtr)[i][j]*s_hat + (*_cScatPtr)[i][j]*s_hat*s_hat);
      this->operator()(i,j)*=(mass*mass-_sAdler0)/_sNorm;
    }
  }
}

template void KMatrixKPiSFocus::evalMatrixTemplate(const double mass, Spin OrbMom);
template void KMatrixKPiSFocus::evalMatrixTemplate(const complex<double> mass, Spin OrbMom);

void KMatrixKPiSFocus::init2IsoSpin1(){

  (*_aScatPtr)[0][0]=0.79299;
  (*_aScatPtr)[0][1]=0.15040;
  (*_aScatPtr)[1][1]=0.17054;
  (*_aScatPtr)[1][0]=(*_aScatPtr)[0][1];

  (*_bScatPtr)[0][0]=-0.15099;
  (*_bScatPtr)[0][1]=-0.038266;
  (*_bScatPtr)[1][1]=-0.0219;
  (*_bScatPtr)[1][0]=(*_bScatPtr)[0][1];

  (*_cScatPtr)[0][0]=0.00811;
  (*_cScatPtr)[0][1]=0.0022596;
  (*_cScatPtr)[1][1]=0.00085655;
  (*_cScatPtr)[1][0]=(*_cScatPtr)[0][1];

  std::map<int, std::vector<double> > gFactorsMap;
  gFactorsMap[0].push_back(0.31072); //Kpi
  gFactorsMap[0].push_back(-0.02323); //Ketap

  std::vector<double> poleMasses;
  poleMasses.resize(1);
  poleMasses[0]=sqrt(1.7919);

  std::map<int, std::vector<double> >::iterator itgFac;
  for (itgFac=gFactorsMap.begin(); itgFac!=gFactorsMap.end(); ++itgFac){
    std::shared_ptr<KPole> currentPole(new KPole(itgFac->second, poleMasses[itgFac->first]));
    _KPoles.push_back(currentPole);
  }

}

void KMatrixKPiSFocus::init2IsoSpin3(){

  (*_aScatPtr)[0][0]=-0.22147;
  (*_aScatPtr)[0][1]=0.;
  (*_aScatPtr)[1][1]=0.;
  (*_aScatPtr)[1][0]=0.;

  (*_bScatPtr)[0][0]=0.026637;
  (*_bScatPtr)[0][1]=0.;
  (*_bScatPtr)[1][1]=0.;
  (*_bScatPtr)[1][0]=0.;

  (*_cScatPtr)[0][0]=-0.00092057;
  (*_cScatPtr)[0][1]=0.;
  (*_cScatPtr)[1][1]=0.;
  (*_cScatPtr)[1][0]=0.;

  std::map<int, std::vector<double> > gFactorsMap;
  gFactorsMap[0].push_back(0.); //Kpi
  gFactorsMap[0].push_back(0.); //Ketap

  std::vector<double> poleMasses;
   poleMasses.resize(1);
   poleMasses[0]=sqrt(100.);

  std::map<int, std::vector<double> >::iterator itgFac;
  for (itgFac=gFactorsMap.begin(); itgFac!=gFactorsMap.end(); ++itgFac){
    std::shared_ptr<KPole> currentPole(new KPole(itgFac->second, poleMasses[itgFac->first]));
    _KPoles.push_back(currentPole);
  }

}
