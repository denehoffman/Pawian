//************************************************************************//
//									                                      //
//  Copyright 2020 Meike Kuessner (mkuessner@ep1.rub.de)		          //
//  	      	   - Ruhr-Universität Bochum 				              //
//									                                      //
//  This file is part of Pawian.					                      //
//									                                      //
//  Pawian is free software: you can redistribute it and/or modify	      //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	              //
//									                                      //
//  Pawian is distributed in the hope that it will be useful,		      //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	      //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	      //
//  GNU General Public License for more details.	      		          //
//									                                      //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	      //
//                                  									  //
//************************************************************************//

#include "ErrLogger/ErrLogger.hh"
#include "PwaDynamics/PhaseSpaceIsobarLUT1D.hh"
#include "Utils/PawianConstants.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "qft++Extension/PawianUtils.hh"

PhaseSpaceIsobarLUT1D::PhaseSpaceIsobarLUT1D(double mass1, double mass2, string type):
  PhaseSpaceIsobar(mass1, mass2)
{
  _name="LUT";
  m_nRe = 0;
  m_nIm = 0;
  m_sLow = 0.0;
  m_sHigh = 0.0;
  m_lutfilepath = type.substr(5, type.length());
  m_sortedByReS = true;
  std::cout << "LUTFile: " << m_lutfilepath << std::endl;
  if (FILE *file = fopen(m_lutfilepath.c_str(), "r")) {
    fclose(file);
  } else {
    Alert << "file: " << m_lutfilepath << " does not exist!!!!" << endmsg;
    exit(1);  
  } 
  loadParams();
}

PhaseSpaceIsobarLUT1D::~PhaseSpaceIsobarLUT1D(){
  m_lutfile.close();//
}

complex<double> PhaseSpaceIsobarLUT1D::factor(const double mass){
  double SqrtSCompl = mass;
  complex<double> result = getFactor(SqrtSCompl);
  complex<double> rho(-result.imag(), 0.); 
  return rho;
}

complex<double> PhaseSpaceIsobarLUT1D::factor(const complex<double> mass){
    Alert << "Provided LUT is one dimensional in s. 2D is work in progress. Use LUT instead of LUTRe for now: " << m_lutfilepath << endmsg;
    return complex<double>(0., 0.);
    exit(1);
}

complex<double> PhaseSpaceIsobarLUT1D::breakUpMom(const double mass){
  complex<double> currentQ=-2.*ChewM(mass)/mass;
  complex<double> result(currentQ.imag(), 0.);
  return result;  
}

complex<double> PhaseSpaceIsobarLUT1D::breakUpMom(const complex<double> mass){
    Alert << "Provided LUT is one dimensional in s. 2D is work in progress. Use LUT instead of LUTRe for now: " << m_lutfilepath << endmsg;
    return complex<double>(0., 0.);
    exit(1);
}

complex<double> PhaseSpaceIsobarLUT1D::ChewM(const complex<double> mass){
    Alert << "Provided LUT is one dimensional in s. 2D is work in progress. Use LUT instead of LUTRe for now: " << m_lutfilepath << endmsg;
    return complex<double>(0., 0.);
    exit(1);
}

complex<double> PhaseSpaceIsobarLUT1D::ChewM(const double mass){
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

void PhaseSpaceIsobarLUT1D::CorrectCMForChosenSign(complex<double>& breakUpMom, complex<double>& toChange){
  if((_bumImPartSign > 0 && breakUpMom.imag() < 0) ||
     (_bumImPartSign < 0 && breakUpMom.imag() > 0)){
     toChange = conj(toChange);
  }
}

void PhaseSpaceIsobarLUT1D::cacheFactors(const double mass){
  int massInt100keV=mass*10000.;
  std::map<int, complex<double> >::const_iterator it = _CMCache.find(massInt100keV);
  if( it == _CMCache.end()){
      double SqrtSCompl = mass;
      complex<double> currentCM=getFactor(SqrtSCompl);
      _CMCache[massInt100keV]=currentCM;
    }
}

LUTLine1D PhaseSpaceIsobarLUT1D::readLine(int _lineOffset){
  LUTLine1D result;
  m_lutfile.seekg(0, m_lutfile.beg);
  m_lutfile.seekg(_lineOffset*4*m_dSize);
  char *memblock = new char [m_dSize];
  m_lutfile.read(memblock, m_dSize);
  double *val=reinterpret_cast<double*>(memblock);
  double ReS = *val;
  
  m_lutfile.read(memblock, m_dSize);
  val=reinterpret_cast<double*>(memblock);
  double dummy = *val;


  m_lutfile.read(memblock, m_dSize);
  val=reinterpret_cast<double*>(memblock);
  double Re = *val;
  m_lutfile.read(memblock, m_dSize);
  val=reinterpret_cast<double*>(memblock);
  double Im = *val;
  double req = ReS;
  complex<double> res = complex<double>(Re, Im);
  result.requested = req;
  result.cm = res;
  return result;
}

complex<double> PhaseSpaceIsobarLUT1D::interpolateCM(LUTLine1D _P1, LUTLine1D _P2, double _s){
  complex<double> result;

  double resultRe =  (_P1.cm.real()*(_P2.requested-_s) + _P2.cm.real()*(_s-_P1.requested) ) / (_P2.requested - _P1.requested);
  double resultIm =  (_P1.cm.imag()*(_P2.requested-_s) + _P2.cm.imag()*(_s-_P1.requested) ) / (_P2.requested - _P1.requested);

  result = complex<double>(resultRe, resultIm);

  return result;
}

complex<double> PhaseSpaceIsobarLUT1D::getFactor(double _s){
  
  if(_s < m_sLow || _s > m_sHigh){
      Alert << "LUT1D: requested CM out of boundaries in file " << m_lutfilepath << endmsg;
      Alert << "m_sLow " << m_sLow << " _s: " << _s << " m_sHigh " << m_sHigh << endmsg;
      return complex<double>(0., 0.);
      exit(1);
  }
  
  LUTLine1D P1;
  LUTLine1D P2;
  int lineOffset = 0;
  double deltaRe = (m_sHigh-m_sLow)/m_nRe;
  
  int nRe = round((_s-m_sLow)/deltaRe);
  if(nRe == 0 || nRe == m_nRe){
    if(abs(_s) > 1.e-15){
      Alert << "Range error while interpolating" << endmsg;
      exit(1);
    }
  }
  
  lineOffset = nRe+1;

  LUTLine1D PointP = readLine(lineOffset);
  if((PointP.requested <= _s)){
    P1 = PointP;
    P2 = readLine(lineOffset+1);
  }
  else if((PointP.requested > _s)){
    P2 = PointP;
    P1 = readLine(lineOffset-1);
  }
  else{
    Alert << "Unforeseen error while interpolating" << endmsg;
    exit(1);
  }
  complex<double> result = interpolateCM(P1, P2, _s);
  
  return result;       
}


void PhaseSpaceIsobarLUT1D::loadParams(){
  m_lutfile.open(m_lutfilepath.c_str(), std::ios::binary | std::ifstream::in);
  if (!m_lutfile.is_open()){
    Alert << "LUT: file not readable: " << m_lutfilepath << endmsg;
    exit(1);
  }
  m_lutfile.seekg(0, m_lutfile.end);
  long length = m_lutfile.tellg();
  m_lutfile.seekg(0, m_lutfile.beg);
  double firstRe = 0.0;
  double lastRe = 0.0;
  int dSize = sizeof(lastRe);
  m_dSize = dSize;
  long lines = length/(4*dSize);
  char *memblock = new char [dSize];
  m_lutfile.read(memblock, dSize);
  double *res=reinterpret_cast<double*>(memblock);
  firstRe = *res;
  m_lutfile.seekg(length-4*dSize);
  m_lutfile.read(memblock, dSize);
  res=reinterpret_cast<double*>(memblock);
  lastRe = *res;
  m_sLow = firstRe;
  m_sHigh = lastRe;
  std::cout << "DEBUG-- sLow: " << m_sLow << " sHigh: " << m_sHigh << " nLines: " << lines << std::endl;
  m_nRe = lines;
  
  return;
}

