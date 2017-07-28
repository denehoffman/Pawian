//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universit??t Bochum 				  //
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
#include <dlfcn.h>

extern "C" {
  void (*jl_initPWA)() =NULL;
  void (*jl_atexit_hookPWA)(int) =NULL;
  jl_value_t* (*jl_box_float64PWA)(double)=NULL;
  jl_module_t* (*jl_loadPWA)(const char*)=NULL;
  jl_function_t* (*jl_get_functionPWA)(jl_module_t*, std::string)=NULL;
  jl_array_t* (*jl_callPWA)(jl_function_t*, jl_value_t**, int)=NULL;
  double* (*jl_array_dataPWA)(jl_array_t*)=NULL;
  bool (*jl_exception_occurredPWA)(void)=NULL;
}

PhaseSpaceIsobarBBUnstable::PhaseSpaceIsobarBBUnstable(double mass1, double mass2, std::string type):
  PhaseSpaceIsobar(mass1, mass2)
{
  obj = dlopen ("libjulia.so", RTLD_LAZY);
  if (!obj) {
    std::string err = dlerror();
    fprintf(stderr, "Fatal in PhaseSpaceIsobarBBUnstable::PhaseSpaceIsobarBBUnstable: cannot load library %s\n", err.c_str());
         exit(1);
      }

  m_pathToModule+= getenv("TOP_DIR");
  m_pathToModule+="/PwaDynamics/basdevant-berger-cm.jl";

  jl_initPWA = (void (*)()) dlsym (obj, "jl_init");
  jl_atexit_hookPWA  = (void (*)(int)) dlsym (obj, "jl_atexit_hook");
  jl_box_float64PWA =(jl_value_t* (*) (double) ) dlsym (obj, "jl_box_float64");
  jl_loadPWA = (jl_module_t* (*)(const char*)) dlsym (obj, "jl_load");
  jl_get_functionPWA = (jl_function_t* (*)(jl_module_t*, std::string)) dlsym (obj, "jl_get_function");
  jl_callPWA = (jl_array_t* (*)(jl_function_t*, jl_value_t**, int)) dlsym (obj, "jl_call");
  jl_exception_occurredPWA = (bool (*)(void)) dlsym (obj, "jl_exception_occurred");
  jl_array_dataPWA = (double* (*) (jl_array_t*)) dlsym (obj, "jl_array_data");
  jl_initPWA();
  //  jl_init();

  m_srealHO = jl_box_float64PWA(0.0);
  m_simagHO = jl_box_float64PWA(0.0);
  m_m1_1HO = jl_box_float64PWA(0.0);
  m_m1_2HO = jl_box_float64PWA(0.0);
  m_m2_1HO = jl_box_float64PWA(0.0);
  m_m2_2HO = jl_box_float64PWA(0.0);
  m_mR1HO = jl_box_float64PWA(0.0);
  m_mR2HO = jl_box_float64PWA(0.0);
  m_f1HO = jl_box_float64PWA(0.0);
  m_f2HO = jl_box_float64PWA(0.0);
  m_epsilonHO = jl_box_float64PWA(0.0);
  loadModule(m_pathToModule.c_str());
  double RhoPiArgs[9] = {0.14, 0.14, 0.14, 0.14, 0.758287544, 0.442718872, 0.14, 0.0, 0.0};
  setArgs(RhoPiArgs);
  //fillMap();
}

PhaseSpaceIsobarBBUnstable::~PhaseSpaceIsobarBBUnstable(){
   jl_atexit_hookPWA(0);
}

complex<double> PhaseSpaceIsobarBBUnstable::factor(const double mass){
  complex<double> massSqrCompl(mass*mass, 1.e-10); // for real s: expansion to s=0 from 1st quadrant
  double sHO[] = {std::real(massSqrCompl), std::imag(massSqrCompl)};
  double* res = computeFactor(sHO);
  complex<double> result(res[0], res[1]);
  return imag(result);
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
  complex<double> resC (res[0], res[1]); 
  complex<double> rho = imag(resC); 
  return rho;
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
  complex<double> massSqrCompl(mass*mass, 1.e-10); // for real s: expansion to s=0 from 1st quadrant
  double sHO[] = {real(massSqrCompl), imag(massSqrCompl)};
  double* res = computeFactor(sHO);
  complex<double> result(res[0], res[1]);
  return result;  
}

complex<double> PhaseSpaceIsobarBBUnstable::ChewM(const complex<double> mass){
  complex<double> s=mass*mass;
  double sHO[] = {real(s), imag(s)};
  double* res = computeFactor(sHO);
  complex<double> result (res[0], res[1]);
  complex<double> momReid = imag(result)*mass/2.0;
  CorrectForChosenSign(momReid, result);
  return result;
}

bool PhaseSpaceIsobarBBUnstable::loadModule(std::string _moduleDefinitionPath){
  bool result = false;
  m_mod = (jl_module_t*) jl_loadPWA(_moduleDefinitionPath.c_str());
  m_func = jl_get_functionPWA(m_mod, "bbcm");
  if (!jl_exception_occurredPWA())result = true;
  return result;
}

double* PhaseSpaceIsobarBBUnstable::computeFactor(double* _inS){
  m_srealHO = jl_box_float64PWA(_inS[0]);
  m_simagHO = jl_box_float64PWA(_inS[1]);
  jl_value_t* allargs[] = {m_srealHO, m_simagHO, m_m1_1HO, m_m1_2HO, m_m2_1HO, m_m2_2HO, m_mR1HO, m_f1HO, m_mR2HO, m_f2HO, m_epsilonHO};
  jl_array_t *ret = (jl_array_t*) jl_callPWA(m_func, allargs, 11);
  double *retVal = (double*) jl_array_dataPWA(ret);
  return retVal;
}

double* PhaseSpaceIsobarBBUnstable::computeFactor(double _inSreal){
  m_srealHO = jl_box_float64PWA(_inSreal);
  m_simagHO = jl_box_float64PWA(0.0);
  jl_value_t* allargs[] = {m_srealHO, m_simagHO, m_m1_1HO, m_m1_2HO, m_m2_1HO, m_m2_2HO, m_mR1HO, m_f1HO, m_mR2HO, m_f2HO, m_epsilonHO};
  jl_array_t *ret = (jl_array_t*) jl_callPWA(m_func, allargs, 11);
  double *retVal = (double*) jl_array_dataPWA(ret);
  return retVal;
}

double* PhaseSpaceIsobarBBUnstable::computeFactor(){
  jl_value_t* allargs[] = {m_srealHO, m_simagHO, m_m1_1HO, m_m1_2HO, m_m2_1HO, m_m2_2HO, m_mR1HO, m_f1HO, m_mR2HO, m_f2HO, m_epsilonHO};
  jl_array_t *ret = (jl_array_t*) jl_callPWA(m_func, allargs, 11);
 double *retVal = (double*) jl_array_dataPWA(ret);
  return retVal;
}

bool PhaseSpaceIsobarBBUnstable::setArgs(double* _inArgs){
  bool result = false;
  m_m1_1HO = jl_box_float64PWA(_inArgs[0]);
  m_m1_2HO = jl_box_float64PWA(_inArgs[1]);
  m_m2_1HO = jl_box_float64PWA(_inArgs[2]);
  m_m2_2HO = jl_box_float64PWA(_inArgs[3]);
  m_mR1HO = jl_box_float64PWA(_inArgs[4]);
  m_f1HO = jl_box_float64PWA(_inArgs[5]);
  m_mR2HO = jl_box_float64PWA(_inArgs[6]);
  m_f2HO = jl_box_float64PWA(_inArgs[7]);
  m_epsilonHO = jl_box_float64PWA(_inArgs[8]); 
  if (!jl_exception_occurredPWA())result = true;
  return result;
}

bool PhaseSpaceIsobarBBUnstable::setS(double* _inS){
  bool result = false;
  m_srealHO = jl_box_float64PWA(_inS[0]);
  m_simagHO = jl_box_float64PWA(_inS[1]);
  if (!jl_exception_occurredPWA())result = true;
  return result;
}



