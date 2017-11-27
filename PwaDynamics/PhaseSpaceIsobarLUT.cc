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
#include "PwaDynamics/PhaseSpaceIsobarLUT.hh"
#include "Utils/PawianConstants.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++Extension/PawianUtils.hh"

PhaseSpaceIsobarLUT::PhaseSpaceIsobarLUT(double mass1, double mass2, string type):
  PhaseSpaceIsobar(mass1, mass2)
{
  m_nRe = 0;
  m_nIm = 0;
  m_sLow = complex<double>(0.0, 0.0);
  m_sHigh = complex<double>(0.0, 0.0);
  m_lutfilepath = type.substr(3, type.length());
  m_sortedByReS = true;
  std::cout << "LUTFile: " << m_lutfilepath << std::endl;
  
  loadParams();
}

PhaseSpaceIsobarLUT::~PhaseSpaceIsobarLUT(){
  //
}

complex<double> PhaseSpaceIsobarLUT::factor(const double mass){
  complex<double> SqrtSCompl(mass, 0.0);
  complex<double> result = getFactor(SqrtSCompl);
  complex<double> rho(result.imag(), 0.); 
  return rho;
}

complex<double> PhaseSpaceIsobarLUT::breakUpMom(const double mass){
  complex<double> SqrtSCompl(mass, 0.0);
  complex<double> result = getFactor(SqrtSCompl);
  complex<double> momReid(imag(result)*mass/2.0, 0.);
  return momReid;
}

complex<double> PhaseSpaceIsobarLUT::factor(const complex<double> mass){
  complex<double> result = getFactor(mass);
  complex<double> rho(result.imag(), 0.); 
  return rho;
}

complex<double> PhaseSpaceIsobarLUT::breakUpMom(const complex<double> mass){
  complex<double> result = getFactor(mass);
  complex<double> momReid= imag(result)*mass/2.0;
  return momReid;
}

complex<double> PhaseSpaceIsobarLUT::ChewM(const double mass){
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

complex<double> PhaseSpaceIsobarLUT::ChewM(const complex<double> mass){
  complex<double> result = getFactor(mass);
  complex<double> momReid = result.imag()*mass/2.;
  CorrectCMForChosenSign(momReid, result); 
  return result;
}

void PhaseSpaceIsobarLUT::CorrectCMForChosenSign(complex<double>& breakUpMom, complex<double>& toChange){
  if((_bumImPartSign > 0 && breakUpMom.imag() < 0) ||
     (_bumImPartSign < 0 && breakUpMom.imag() > 0)){
     toChange = conj(toChange);
  }
}

void PhaseSpaceIsobarLUT::cacheFactors(const double mass){
  int massInt100keV=mass*10000.;
  std::map<int, complex<double> >::const_iterator it = _CMCache.find(massInt100keV);
  if( it == _CMCache.end()){
      complex<double> SqrtSCompl(mass, 0.0);
      complex<double> currentCM=getFactor(SqrtSCompl);
      _CMCache[massInt100keV]=currentCM;
    }
}


complex<double> PhaseSpaceIsobarLUT::getFactor(complex<double> _s){
  if(_s.real() < m_sLow.real() || _s.real() > m_sHigh.real() || _s.imag() < m_sHigh.imag() || _s.imag() > m_sLow.imag()){
   return complex<double>(0., 0.);
  }
  int lineOffset = 0;
  double deltaRe = (m_sHigh.real()-m_sLow.real())/m_nRe;
  double deltaIm = (m_sHigh.imag()-m_sLow.imag())/m_nIm;
  int dSize = sizeof(deltaRe);
  int nRe = (_s.real()-m_sLow.real())/deltaRe;
  int nIm = (_s.imag()-m_sLow.imag())/deltaIm;
  if(m_nIm == 1) nIm = 0;
  lineOffset = m_nIm*nRe+nIm;
  m_lutfile.open(m_lutfilepath.c_str(), std::ios::binary | std::ifstream::in);
  m_lutfile.seekg(0, m_lutfile.beg);
  m_lutfile.seekg(lineOffset*4*dSize);
  char *memblock = new char [dSize];
  m_lutfile.read(memblock, dSize);
  double *res=reinterpret_cast<double*>(memblock);
  double ReS = *res;
  m_lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  double ImS = *res;
  m_lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  double Re = *res;
  m_lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  double Im = *res;
  m_lutfile.close();
  complex<double> reqested = complex<double>(ReS, ImS);
  complex<double> result = complex<double>(Re, Im);

  if(_s.real() > ReS && nRe < m_nRe) nRe++;
  if(_s.real() < ReS && nRe > 0) nRe--;
  if(_s.imag() > ImS && nIm > 0) nIm--;
  if(_s.imag() < ImS && nIm < m_nIm) nIm++;
  lineOffset = m_nIm*nRe+nIm;
  m_lutfile.open(m_lutfilepath.c_str(), std::ios::binary | std::ifstream::in);
  m_lutfile.seekg(0, m_lutfile.beg);
  m_lutfile.seekg(lineOffset*4*dSize);
  m_lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  double ReS2 = *res;
  m_lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  double ImS2 = *res;
  m_lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  double Re2 = *res;
  m_lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  double Im2 = *res;
  m_lutfile.close();
  complex<double> result2 = complex<double>(Re2, Im2);
  double resRe;
  if(ReS!=ReS2) resRe = Re + (_s.real() - ReS) * ((Re2 - Re)/(ReS2 - ReS));
  else resRe = Re2;
  double resIm;
  if(ImS!=ImS2) resIm = Im + (_s.imag() - ImS) * ((Im2 - Im)/(ImS2 - ImS));
  else resIm = Im2;
  result = complex<double>(resRe, resIm);
  //std::cout << "Req: " << _s << " Result: " << result << std::endl;
  return result;
}


void PhaseSpaceIsobarLUT::loadParams(){
  m_lutfile.open(m_lutfilepath.c_str(), std::ios::binary | std::ifstream::in);
  m_lutfile.seekg(0, m_lutfile.beg);
  m_lutfile.seekg(0, m_lutfile.end);
  long size = m_lutfile.tellg();
  m_lutfile.seekg(0, m_lutfile.beg);
  double firstRe = 0.0;
  double lastRe = 0.0;
  double firstIm = 0.0;
  double lastIm = 0.0;
  int dSize = sizeof(lastRe);
  long lines = size/(4*dSize);
  char *memblock = new char [dSize];
  m_lutfile.read(memblock, dSize);
  double *res=reinterpret_cast<double*>(memblock);
  firstRe = *res;
  m_lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  firstIm = *res;
  m_lutfile.seekg(size-4*dSize);
  m_lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  lastRe = *res;
  m_lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  lastIm = *res;
  m_sLow = complex<double>(firstRe, firstIm);
  m_sHigh = complex<double>(lastRe, lastIm);
  //std::cout << "DEBUG-- sLow: " << m_sLow << " sHigh: " << m_sHigh << std::endl;
  m_lutfile.seekg(0, m_lutfile.beg);
  int offset = 1;
  int nReFoundInBlock = 1;
  bool nReBlockFound = false;
  while((!nReBlockFound) && offset < size/(4*dSize)){
    m_lutfile.seekg(offset*4*dSize);
    m_lutfile.read(memblock, dSize);
    res=reinterpret_cast<double*>(memblock);
    if(abs(*res - firstRe) > 1.e-15) nReBlockFound = true;
    else nReFoundInBlock++;
    offset++;
  }
  m_nIm = nReFoundInBlock;
  m_nRe = lines/m_nIm;
  m_lutfile.close();
  return;
}

