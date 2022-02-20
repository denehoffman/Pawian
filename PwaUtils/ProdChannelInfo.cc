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

// ProdChannelInfo class definition file. -*- C++ -*-
// Copyright 2015 Bertram Kopf

#include "PwaUtils/ProdChannelInfo.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "ErrLogger/ErrLogger.hh"

ProdChannelInfo::ProdChannelInfo() :
  _isProdChannel(false)
  ,_isFormation(false)
  , _formationParticle(0)
  , _withProBarrier(false)
  ,_prodBarrierType("non")
  ,_formationDynType("non")
  ,_qRProd(0.1973)
  ,_m0decRadM1Production(1.5)
{
}

ProdChannelInfo::ProdChannelInfo(std::string& stringFromParser) :
  _isProdChannel(true)
  ,_isFormation(false)
  ,_formationParticle(0)
  , _withProBarrier(false)
  ,_prodBarrierType("non")
  ,_qRProd(0.1973)
  ,_m0decRadM1Production(1.5)
{
  std::stringstream strStrFromParser;
  strStrFromParser << stringFromParser;
  
  std::string firstParticleStr;
  std::string secondParticleStr;
  if(! (strStrFromParser >> firstParticleStr) || !(strStrFromParser >> secondParticleStr) ){
    Alert <<"at least two arguments for the production string are required, namely particleName1 and particleName2" << endmsg;
    exit(1);
  }
  
  Particle* firstParticle = GlobalEnv::instance()->particleTable()->particle(firstParticleStr);
  if( 0==firstParticle){
    Alert << "particle\t" << firstParticleStr << "\tdoes not exist in pdtTable" << endmsg;
    exit(1);
  }

  Particle* secondParticle=0;
  if(secondParticleStr=="formation"){
    _isFormation=true;
    _formationParticle=firstParticle;
    //    secondParticle=firstParticle;
    secondParticle=GlobalEnv::instance()->particleTable()->particle("photon"); //just needed to make the compiler happy; must be changed
  }
  else{
    secondParticle = GlobalEnv::instance()->particleTable()->particle(secondParticleStr);
    if( 0==secondParticle){
      Alert << "particle\t" << secondParticleStr << "\tdoes not exist in pdtTable" << endmsg;
      exit(1);
    }
  }

  _prodPair=make_pair(firstParticle, secondParticle);
  
  //look for additional options
  std::string typeStr;
  if( strStrFromParser >> typeStr){
    _withProBarrier=true;
    if(typeStr=="BlattWBarrier" || typeStr=="BlattWBarrierTensor"){
      _prodBarrierType=typeStr;
      std::string qRStr;
      if( strStrFromParser >> qRStr){
	_qRProd=atof(qRStr.c_str());
      }
    }
    else if(typeStr=="RadM1"){
      _prodBarrierType="RadM1";
      std::string m0DecStr;
      if( strStrFromParser >> m0DecStr){
      _m0decRadM1Production=atof(m0DecStr.c_str());
      }
    }
    else if(typeStr=="RadM1KEDR"){
      _prodBarrierType="RadM1KEDR";
      std::string m0DecStr;
      if( strStrFromParser >> m0DecStr){
      _m0decRadM1Production=atof(m0DecStr.c_str());
      }
    }
    else if(typeStr=="FormPol2"){
      if(!_isFormation){
	Alert << "dynamics " << typeStr
	    << "\tis only allowed for the formation process" << endmsg;
      }
       _formationDynType=typeStr;
    }
    else if(typeStr=="SExpDynamics"){
      if(!_isFormation){
        Alert << "dynamics " << typeStr
            << "\tis only allowed for the formation process" << endmsg;
      }
    _prodBarrierType=typeStr;
    _formationDynType=typeStr;
    }
    else{
      Alert << "production barrier type: " << typeStr << " doesn't exist for PAWIAN!!!" 
	    << "\nonly BlattWBarrier, BlattWBarrierTensor, RadM1 and RadM1KEDR are supported" << endmsg;
      exit(1);
    }
  }

}


ProdChannelInfo::~ProdChannelInfo(){
}


void ProdChannelInfo::print(std::ostream& os) const{
  os << "\n***ProdChannelInfo***"; 
  if(!_isProdChannel){
    os << "\nIt's not a production channel!!!" << endl;
  }

  os << "\nparticle1: " << _prodPair.first->name() << "\tparticle2: " << _prodPair.second->name()
     << "\nwith production barrier: " << _withProBarrier << std::endl;

  if (_withProBarrier){
    os << "prod barrier type: " << _prodBarrierType << std::endl;
    if(_prodBarrierType=="BlattWBarrier" || _prodBarrierType=="BlattWBarrierTensor") os << "with qR = "  << _qRProd << std::endl;
    else if(_prodBarrierType=="RadM1" || _prodBarrierType=="RadM1KEDR") os << "with m0decRadM1Production = "  << _m0decRadM1Production << std::endl;
  }

}

