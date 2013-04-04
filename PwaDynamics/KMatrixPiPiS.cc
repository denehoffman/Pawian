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

#include "PwaDynamics/KMatrixPiPiS.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceIsobar.hh"
#include "PwaDynamics/PhaseSpace4Pi.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++/matrix/IdentityMatrix.hh"
#include <boost/multi_array.hpp>


KMatrixPiPiS::KMatrixPiPiS() :
  KMatrixSlowAdlerCorRel(5)
{
  initASParam1900();
}


KMatrixPiPiS::~KMatrixPiPiS(){
}

void  KMatrixPiPiS::initASParam1900(){
  const double piMass=0.1349766;
  const double KplusMass=0.493677;
  const double K0Mass= 0.497614;
  const double etaMass= 0.547853;
  const double etaprimeMass= 0.95778;

  std::map<int, std::vector<double> > gFactorsMap;
  std::vector<double> poleMasses;
  poleMasses.resize(5);
  poleMasses[0]=0.651;
  poleMasses[1]=1.2036;
  poleMasses[2]=1.55817;
  poleMasses[3]=1.21;
  poleMasses[4]=1.82206;

  gFactorsMap[0].push_back(0.22889); //pipi
  gFactorsMap[0].push_back(-0.55377); //KK
  gFactorsMap[0].push_back(0.); //4pi
  gFactorsMap[0].push_back(-0.39899); //eta eta
  gFactorsMap[0].push_back(-0.34639); //eta eta'

  gFactorsMap[1].push_back(0.94128); //pipi
  gFactorsMap[1].push_back(0.55095); //KK
  gFactorsMap[1].push_back(0.); //4pi
  gFactorsMap[1].push_back(0.39065); //eta eta
  gFactorsMap[1].push_back(0.31503); //eta eta'

  gFactorsMap[2].push_back(0.36856); //pipi
  gFactorsMap[2].push_back(0.23888); //KK
  gFactorsMap[2].push_back(0.55639); //4pi
  gFactorsMap[2].push_back(0.18340); //eta eta
  gFactorsMap[2].push_back(0.18681); //eta eta'

  gFactorsMap[3].push_back(0.3365); //pipi
  gFactorsMap[3].push_back(0.40907); //KK
  gFactorsMap[3].push_back(0.85679); 
  gFactorsMap[3].push_back(0.19906); //eta eta
  gFactorsMap[3].push_back(-0.00984); //eta eta'

  gFactorsMap[4].push_back(0.18171); //pipi
  gFactorsMap[4].push_back(-0.17558); //KK
  gFactorsMap[4].push_back(-0.79658); //4pi
  gFactorsMap[4].push_back(-0.00355); //eta eta
  gFactorsMap[4].push_back(0.22358); //eta eta'  

  boost::shared_ptr<AbsPhaseSpace> pipiPhp(new PhaseSpaceIsobar(piMass, piMass));
  boost::shared_ptr<AbsPhaseSpace> kkPhp(new PhaseSpaceIsobar(KplusMass, K0Mass));
  boost::shared_ptr<AbsPhaseSpace> pipipipiPhp(new PhaseSpace4Pi());
  boost::shared_ptr<AbsPhaseSpace> etaetaPhp(new PhaseSpaceIsobar(etaMass, etaMass));
  boost::shared_ptr<AbsPhaseSpace> etaetapPhp(new PhaseSpaceIsobar(etaMass, etaprimeMass));

  _phpVecs.push_back(pipiPhp);
  _phpVecs.push_back(kkPhp);
  _phpVecs.push_back(pipipipiPhp);
  _phpVecs.push_back(etaetaPhp);
  _phpVecs.push_back(etaetapPhp);

  //  boost::shared_ptr< boost::multi_array<double, 2> > fScatMatr(new boost::multi_array<double, 2>(boost::extents[gFactors.size()][gFactors.size()]));

  _fScatPtr=boost::shared_ptr< boost::multi_array<double, 2> >(new boost::multi_array<double, 2>(boost::extents[NumRows()][NumCols()]));

  for (int i=0; i<NumRows(); ++i){
    for (int j=0; j<NumCols(); ++j){
      (*_fScatPtr)[i][j]=0.;
    }
  }
  (*_fScatPtr)[0][0]=0.23399;
  (*_fScatPtr)[0][1]=0.15044;
  (*_fScatPtr)[0][2]=-0.20545;
  (*_fScatPtr)[0][3]=0.32825;
  (*_fScatPtr)[0][4]=0.35412;

  (*_fScatPtr)[1][0]=(*_fScatPtr)[0][1];
  (*_fScatPtr)[2][0]=(*_fScatPtr)[0][2];
  (*_fScatPtr)[3][0]=(*_fScatPtr)[0][3];
  (*_fScatPtr)[4][0]=(*_fScatPtr)[0][4];

  std::vector< boost::shared_ptr<KPole> > kPoles;

  std::map<int, std::vector<double> >::iterator itgFac;
  for (itgFac=gFactorsMap.begin(); itgFac!=gFactorsMap.end(); ++itgFac){
    boost::shared_ptr<KPole> currentPole(new KPole(itgFac->second, poleMasses[itgFac->first]));
    _KPoles.push_back(currentPole);
  }

  _s0Scat=-3.92637;
  _sAdler0=-0.15;
  _sAdler=1.;


}
