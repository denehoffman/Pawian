//************************************************************************//
//									  //
//  Copyright 2017 Markus Kuhlmann (mkuhlmann@ep1.rub.de)		  //
//  	      	   - Ruhr-Universität Bochum 				  //
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

#include "ErrLogger/ErrLogger.hh"
#include "PwaDynamics/PhaseSpaceIsobarBBUnstable.hh"
#include "Utils/PawianConstants.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++Extension/PawianUtils.hh"

PhaseSpaceIsobarBBUnstable::PhaseSpaceIsobarBBUnstable(double mass1, double mass2, std::string type):
  PhaseSpaceIsobar(mass1, mass2)
{
  m_knownCombinations.clear();
  m_pathToModule+= getenv("TOP_DIR");
  m_pathToModule+="/PwaDynamics/basdevant-berger-cm.jl";
  jl_init();
  m_srealHO = jl_box_float64(0.0);
  m_simagHO = jl_box_float64(0.0);
  m_m1_1HO = jl_box_float64(0.0);
  m_m1_2HO = jl_box_float64(0.0);
  m_m2_1HO = jl_box_float64(0.0);
  m_m2_2HO = jl_box_float64(0.0);
  m_mR1HO = jl_box_float64(0.0);
  m_mR2HO = jl_box_float64(0.0);
  m_f1HO = jl_box_float64(0.0);
  m_f2HO = jl_box_float64(0.0);
  m_epsilonHO = jl_box_float64(0.0);
  loadModule(m_pathToModule.c_str());
  fillMap();
  std::map<std::string, double*>::iterator it = m_knownCombinations.find(type);
  if(it != m_knownCombinations.end()){
    setArgs(m_knownCombinations.find(type)->second);
  }
  else{
    double paramsGiven[9];
    std::stringstream paramStream(type.substr(11));
    for(unsigned int i = 0; i < 9; i++){
      paramStream >> paramsGiven[i];
    }
    setArgs(paramsGiven);
  }
  
  
}

PhaseSpaceIsobarBBUnstable::~PhaseSpaceIsobarBBUnstable(){
  jl_atexit_hook(0);
}

complex<double> PhaseSpaceIsobarBBUnstable::factor(const double mass){
  complex<double> massSqrCompl(mass*mass, 1.e-10);
  double sHO[] = {std::real(massSqrCompl), std::imag(massSqrCompl)};
  double* res = computeFactor(sHO);
  complex<double> result( res[1], 0);
  return result;
}

complex<double> PhaseSpaceIsobarBBUnstable::breakUpMom(const double mass){
  complex<double> s=mass*mass;
  double sHO[] = {real(s), 1.e-10};
  double* res = computeFactor(sHO);
  complex<double> result (res[0], res[1]);
  complex<double> momReid = imag(result)*mass/2.0;
  CorrectForChosenSign(momReid, momReid);
  return momReid;
}

complex<double> PhaseSpaceIsobarBBUnstable::factor(const complex<double> mass){
  complex<double> s=mass*mass;
  double sHO[] = {real(s), imag(s)};
  double* res = computeFactor(sHO);
  complex<double> resC (res[1], 0.); 
  return resC;
}

complex<double> PhaseSpaceIsobarBBUnstable::breakUpMom(const complex<double> mass){
  complex<double> s=mass*mass;
  double sHO[] = {real(s), imag(s)};
  double* res = computeFactor(sHO);
  complex<double> result (res[0], res[1]);
  complex<double> momReid = imag(result)*mass/2.0;
  CorrectForChosenSign(momReid, momReid);
  return momReid;
}

complex<double> PhaseSpaceIsobarBBUnstable::ChewM(const double mass){
  complex<double> result(0.,0.);
  int massInt100keV=mass*10000.;
  std::map<int, complex<double> >::const_iterator it = _CMCache.find(massInt100keV);
  if(it != _CMCache.end()){
    result=it->second;
  }
  else{
    WarningMsg << "_CMCache not found for mass/100keV: " << massInt100keV 
	    << "\t cach it now!!!" << endmsg;
    //    exit(1);
    cacheFactors(mass);
    result=_CMCache.at(massInt100keV); 
  }

  return result;
}

complex<double> PhaseSpaceIsobarBBUnstable::ChewM(const complex<double> mass){
  complex<double> massSqrCompl=mass*mass;
  double sHO[] = {real(massSqrCompl), imag(massSqrCompl)};
  double* res = computeFactor(sHO);  
  complex<double> result (res[0], res[1]);
  complex<double> momReid = imag(result)*mass/2.0;
  CorrectForChosenSign(momReid, result);
  return result;
}

bool PhaseSpaceIsobarBBUnstable::loadModule(std::string _moduleDefinitionPath){
  bool result = false;
  m_mod = (jl_module_t*) jl_load(_moduleDefinitionPath.c_str());
  m_func = jl_get_function(m_mod, "bbcm");
  if (!jl_exception_occurred())result = true;
  return result;
}

double* PhaseSpaceIsobarBBUnstable::computeFactor(double* _inS){
  m_srealHO = jl_box_float64(_inS[0]);
  m_simagHO = jl_box_float64(_inS[1]);
  jl_value_t* allargs[] = {m_srealHO, m_simagHO, m_m1_1HO, m_m1_2HO, m_m2_1HO, m_m2_2HO, m_mR1HO, m_f1HO, m_mR2HO, m_f2HO, m_epsilonHO};
  jl_array_t *ret = (jl_array_t*) jl_call(m_func, allargs, 11);
  double *retVal = (double*) jl_array_data(ret);
  return retVal;
}

complex<double> PhaseSpaceIsobarBBUnstable::computeFactor(complex<double> inS){
  std::cout << "HERE " << inS << std::endl;
  m_srealHO = jl_box_float64(inS.real());
  m_simagHO = jl_box_float64(inS.imag());
  jl_value_t* allargs[] = {m_srealHO, m_simagHO, m_m1_1HO, m_m1_2HO, m_m2_1HO, m_m2_2HO, m_mR1HO, m_f1HO, m_mR2HO, m_f2HO, m_epsilonHO};
  jl_array_t *ret = (jl_array_t*) jl_call(m_func, allargs, 11);
  double *retVal = (double*) jl_array_data(ret);
  complex<double> result( retVal[0], retVal[1]);
  return result;
}

double* PhaseSpaceIsobarBBUnstable::computeFactor(double _inSreal){
  m_srealHO = jl_box_float64(_inSreal);
  m_simagHO = jl_box_float64(0.0);
  jl_value_t* allargs[] = {m_srealHO, m_simagHO, m_m1_1HO, m_m1_2HO, m_m2_1HO, m_m2_2HO, m_mR1HO, m_f1HO, m_mR2HO, m_f2HO, m_epsilonHO};
  jl_array_t *ret = (jl_array_t*) jl_call(m_func, allargs, 11);
  double *retVal = (double*) jl_array_data(ret);
  return retVal;
}

double* PhaseSpaceIsobarBBUnstable::computeFactor(){
  jl_value_t* allargs[] = {m_srealHO, m_simagHO, m_m1_1HO, m_m1_2HO, m_m2_1HO, m_m2_2HO, m_mR1HO, m_f1HO, m_mR2HO, m_f2HO, m_epsilonHO};
  jl_array_t *ret = (jl_array_t*) jl_call(m_func, allargs, 11);
  double *retVal = (double*) jl_array_data(ret);
  return retVal;
}

bool PhaseSpaceIsobarBBUnstable::setArgs(double* _inArgs){
  bool result = false;
  m_m1_1HO = jl_box_float64(_inArgs[0]);
  m_m1_2HO = jl_box_float64(_inArgs[1]);
  m_m2_1HO = jl_box_float64(_inArgs[2]);
  m_m2_2HO = jl_box_float64(_inArgs[3]);
  m_mR1HO = jl_box_float64(_inArgs[4]);
  m_f1HO = jl_box_float64(_inArgs[5]);
  m_mR2HO = jl_box_float64(_inArgs[6]);
  m_f2HO = jl_box_float64(_inArgs[7]);
  m_epsilonHO = jl_box_float64(_inArgs[8]);
  if (!jl_exception_occurred())result = true;
  return result;
}

bool PhaseSpaceIsobarBBUnstable::setS(double* _inS){
  bool result = false;
  m_srealHO = jl_box_float64(_inS[0]);
  m_simagHO = jl_box_float64(_inS[1]);
  if (!jl_exception_occurred())result = true;
  return result;
}

void PhaseSpaceIsobarBBUnstable::fillMap(){
  
  double RhoPiArgs[9] = {0.13957, 0.13957, 0.13957, 0.13957, 0.758287544, 0.442718872, 0.13957, 0.0, 0.0};
  double RhoKArgs[9] = {0.13957, 0.13957, 0.13957, 0.13957, 0.758287544, 0.442718872, 0.493667, 0.0, 0.0};
  double RhoRhoArgs[9] = {0.13957, 0.13957, 0.13957, 0.13957, 0.758287544, 0.442718872, 0.758287544, 0.442718872, 0.0};
  double KstarKArgs[9] = {0.493667, 0.13957, 0.13957, 0.13957, 0.903880523, 0.425440948, 0.493667, 0.0, 0.0};
  double KstarPiArgs[9] = {0.493667, 0.13957, 0.13957, 0.13957, 0.903880523, 0.425440948, 0.13957, 0.0, 0.0};
  double KstarKstarArgs[9] = {0.493667, 0.13957, 0.493667, 0.13957, 0.903880523, 0.425440948, 0.903880523, 0.425440948, 0.0};
  double KstarRhoArgs[9] = {0.493667, 0.13957, 0.493667, 0.13957, 0.903880523, 0.425440948, 0.758287544, 0.442718872, 0.0};
  m_knownCombinations.insert ( std::pair<std::string,double*>("BBUnstableRhoPi",RhoPiArgs) );
  m_knownCombinations.insert ( std::pair<std::string,double*>("BBUnstableRhoK",RhoKArgs) );
  m_knownCombinations.insert ( std::pair<std::string,double*>("BBUnstableRhoRho",RhoRhoArgs) );
  m_knownCombinations.insert ( std::pair<std::string,double*>("BBUnstableKstarK",KstarKArgs) );
  m_knownCombinations.insert ( std::pair<std::string,double*>("BBUnstableKstarPi",KstarPiArgs) );
  m_knownCombinations.insert ( std::pair<std::string,double*>("BBUnstableKstarKstar",KstarKstarArgs) );
  m_knownCombinations.insert ( std::pair<std::string,double*>("BBUnstableKstarRho",KstarRhoArgs) );
  return;
}

void PhaseSpaceIsobarBBUnstable::cacheFactors(const double mass){
  int massInt100keV=mass*10000.;
  std::map<int, complex<double> >::const_iterator it = _CMCache.find(massInt100keV);
  if( it == _CMCache.end()){
      complex<double> massSqrCompl(mass*mass, 1.e-10);
      complex<double> currentCM=computeFactor(massSqrCompl);
      _CMCache[massInt100keV]=currentCM;
    }
}

