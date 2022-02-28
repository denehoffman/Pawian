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

// centralProdBaseLh class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "ggUtils/centralProdBaseLh.hh"
#include "ggUtils/ggReaction.hh"
#include "ggUtils/GGChannelEnv.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/FsParticleProjections.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"
#include "ConfigParser/ParserBase.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

centralProdBaseLh::centralProdBaseLh(ChannelID channelID) :
   AbsLh(channelID)
   ,_ggChannelEnv(std::static_pointer_cast<GGChannelEnv> (GlobalEnv::instance()->GGChannel(channelID)))
   ,_useProdDynamics(false) {
   initialize();
}


centralProdBaseLh::~centralProdBaseLh()
{;
}


complex<double> centralProdBaseLh::calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, 
						   EvtData* theData){
   complex<double> result(0.,0.);
   return result;
}


complex<double> centralProdBaseLh::calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, 
						   EvtData* theData,
      std::map <std::shared_ptr<const JPCLS>,
      std::vector< std::shared_ptr<AbsXdecAmp> >,
      pawian::Collection::SharedPtrLess > pbarpAmps) {
   complex<double> resultAmp(0.,0.);

   return resultAmp;
}


double centralProdBaseLh::calcEvtIntensity( EvtData* theData, std::shared_ptr<AbsPawianParameters> fitPar){

   double result=0.;

   std::vector< std::shared_ptr<AbsXdecAmp> >::iterator itDecAll;
   for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll){
      (*itDecAll)->calcDynamics(theData);
   }

   const std::vector< std::vector<Spin> >& spinProjections=_fsParticleProjections->spinProjections();

   // BEGIN Outer incoherent sum over spinProjections of FS particles...
   for (unsigned int projId=0; projId<spinProjections.size(); ++projId) {
     for (itDecAll=_decAmps.begin(); itDecAll!=_decAmps.end(); ++itDecAll) {
       //(*itDecAll)->setSpinProjections(projId);
       const std::vector<Spin>& currentSpinProjection=spinProjections.at(projId);
      (*itDecAll)->setSpinProjections(currentSpinProjection, projId);
     }

     std::vector< complex<double> > lamXPomPomAmps;
     std::vector< complex<double> > lamXGamPomAmps;
     
     std::vector<Spin>::iterator itLamX;
     for(itLamX=_lamX.begin(); itLamX!=_lamX.end(); ++itLamX) {
       lamXPomPomAmps.push_back(complex<double>(0.,0.));
       lamXGamPomAmps.push_back(complex<double>(0.,0.));

       std::vector<std::shared_ptr<AbsXdecAmp> >::iterator itDec;
       for( itDec=_pompomAmps.begin(); itDec!=_pompomAmps.end(); ++itDec){
	 complex<double> currentDecAmp=(*itDec)->XdecAmp((*itLamX), theData);
	 lamXPomPomAmps.back() += currentDecAmp;
       }

       for( itDec=_gampomAmps.begin(); itDec!=_gampomAmps.end(); ++itDec){
	 complex<double> currentDecAmp=(*itDec)->XdecAmp((*itLamX), theData);
	 lamXGamPomAmps.back() += currentDecAmp;
       }
       if(_useCohPhasespace) //coherent phases space only for pompom events
	 lamXPomPomAmps.back() += std::polar( fitPar->Value(_CohPhasespaceKey+"Mag"), 
					      fitPar->Value(_CohPhasespaceKey+"Phi") );
      }

      std::vector< complex<double> >::iterator itLamXAmps;
      for(itLamXAmps=lamXPomPomAmps.begin(); itLamXAmps!=lamXPomPomAmps.end(); ++itLamXAmps) {
         result += norm(*itLamXAmps);
      }
      for(itLamXAmps=lamXGamPomAmps.begin(); itLamXAmps!=lamXGamPomAmps.end(); ++itLamXAmps) {
         result += norm(*itLamXAmps);
      }
   } // END Outer incoherent sum over spinProjections of FS particles.

   if(_usePhasespace) result+=fitPar->Value(_phasespaceKey);
   if(_useProdDynamics) result*=fabs(_dyn->eval(theData,0).real());
   result *= fitPar->Value(_channelScaleParam);
   return result;
}

void centralProdBaseLh::print(std::ostream& os) const{

}

void  centralProdBaseLh::initialize(){
   std::vector< std::shared_ptr<AbsDecay> > theDecs = _ggChannelEnv->prodDecayList()->getList();
   std::vector< std::shared_ptr<AbsDecay> >::iterator it;
   int maxJ = 0;
   for (it=theDecs.begin(); it!=theDecs.end(); ++it){
     InfoMsg << "theDecs->name: " << (*it)->name() << endmsg;
     std::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(_channelID, 
										    (*it)->absDecPtr());
     _decAmps.push_back(currentAmp);
      
     bool hasPosCparity=false;
     if ( currentAmp->jpcPtr()->C == -1)  hasPosCparity=true;
     
     
     if(hasPosCparity) _gampomAmps.push_back(currentAmp);
     else _pompomAmps.push_back(currentAmp);
     
     if(currentAmp->jpcPtr()->J > maxJ) maxJ = currentAmp->jpcPtr()->J;
   }
   for(int i = (-1*maxJ); i<= maxJ; i++) _lamX.push_back(Spin(i));
   
   std::string dynString=_ggChannelEnv->parser()->productionDynamics();
   if(dynString=="woFormPol" || dynString=="FormPol0" || dynString=="FormPol1" || dynString=="FormPol2"){
     _useProdDynamics=true;
     std::shared_ptr<FormationDecay> motherFormDec = _ggChannelEnv->reaction()->motherProdDec();
     std::vector<std::string> additionalStringVecDummy;
     motherFormDec->enableDynamics(dynString, additionalStringVecDummy);
     _dyn=motherFormDec->getDynamics();
   }
}

void centralProdBaseLh::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  if(_useProdDynamics) _dyn->fillDefaultParams(fitPar);
  AbsLh::fillDefaultParams(fitPar);
}

void centralProdBaseLh::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  if(_useProdDynamics) _dyn->updateFitParams(fitPar);
  AbsLh::updateFitParams(fitPar);
}
