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

// pbarpReaction class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "pbarpUtils/pbarpReaction.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "pbarpUtils/pbarpStatesLS.hh"

pbarpReaction::pbarpReaction(std::vector<std::pair<Particle*, Particle*> >& prodPairs, int lmax) :
  _lmax(lmax)
{
  boost::shared_ptr<pbarpStatesLS> thepbarpStates(new pbarpStatesLS(lmax));
  std::vector< boost::shared_ptr<const jpcRes> > pbarpJPCStatesAll= thepbarpStates->jpcStates(); 
  std::vector< boost::shared_ptr<const jpcRes> >::const_iterator itJPC;
  std::vector< boost::shared_ptr<const JPCLS> > all_JPCLSs= thepbarpStates->all_JPCLS_States();
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  for(itJPCLS = all_JPCLSs.begin(); itJPCLS != all_JPCLSs.end(); ++itJPCLS){
     boost::shared_ptr<const jpcRes> currentJPC = (*itJPCLS);
     auto minMapIt = _minLMap.find(currentJPC);
     if(minMapIt == _minLMap.end()){
	_minLMap[currentJPC] = (*itJPCLS)->L;
     }
     else{
	if(_minLMap[currentJPC] > (*itJPCLS)->L)
	   _minLMap[currentJPC] = (*itJPCLS)->L;
     }

     _jpcToLMap[currentJPC].push_back((*itJPCLS)->L);
  }

  for(itJPC = pbarpJPCStatesAll.begin(); itJPC!=pbarpJPCStatesAll.end(); ++itJPC){
    bool acceptJPC=false;
    std::vector<std::pair<Particle*, Particle*> >::iterator itPartPairs;
    for (itPartPairs=prodPairs.begin(); itPartPairs!= prodPairs.end(); ++itPartPairs){
      std::string decName=(*itJPC)->name();
      //      std::string decName="";
      boost::shared_ptr<IsobarLSDecay> currentDec(new IsobarLSDecay( (*itJPC),itPartPairs->first, itPartPairs->second, pbarpEnv::instance(), decName));

      if (currentDec->JPCLSAmps().size()>0 &&
	  CheckLDrop((std::string&)(*itPartPairs).first->name(), *itJPC) &&
	  CheckLDrop((std::string&)(*itPartPairs).second->name(), *itJPC)){
	 _prodDecs.push_back(currentDec);

	 //boost::shared_ptr<IsobarTensorDecay> currentTensorDec(new IsobarTensorDecay( (*itJPC),itPartPairs->first, itPartPairs->second, pbarpEnv::instance(), decName));
         //_prodTensorDecs.push_back(currentTensorDec);

	 boost::shared_ptr<IsobarHeliDecay> currentHeliDec(new IsobarHeliDecay( (*itJPC),itPartPairs->first, itPartPairs->second, pbarpEnv::instance(), decName));
	 _prodHeliDecs.push_back(currentHeliDec);
	 acceptJPC=true;
      }
    }
    if(acceptJPC)
       _pbarpJPCs.push_back(*itJPC);
  }

  _pbarpJPCLSs =  extractStates(_pbarpJPCs, all_JPCLSs);

  std::vector< boost::shared_ptr<const JPCLS> > all_pbarpSingletLS = thepbarpStates->singlet_JPCLS_States();
  _pbarpJPCLSsinglet =  extractStates(_pbarpJPCs, all_pbarpSingletLS);
  // fillMap(_pbarpJPCLSsinglet, _prodDecs, _pbarpSingletDecMap);

  std::vector< boost::shared_ptr<const JPCLS> > all_pbarpTriplet0LS = thepbarpStates->triplet0_JPCLS_States();
  _pbarpJPCLStriplet0 =  extractStates(_pbarpJPCs, all_pbarpTriplet0LS);
  // fillMap(_pbarpJPCLStriplet0, _prodDecs, _pbarpTriplet0DecMap);

  std::vector< boost::shared_ptr<const JPCLS> > all_pbarpTripletp1LS = thepbarpStates->tripletp1_JPCLS_States();
  _pbarpJPCLStripletp1 =  extractStates(_pbarpJPCs, all_pbarpTripletp1LS);
  // fillMap(_pbarpJPCLStripletp1, _prodDecs, _pbarpTripletp1DecMap);

  std::vector< boost::shared_ptr<const JPCLS> > all_pbarpTripletm1LS = thepbarpStates->tripletm1_JPCLS_States();
  _pbarpJPCLStripletm1 =  extractStates(_pbarpJPCs, all_pbarpTripletm1LS);
  // fillMap(_pbarpJPCLStripletm1, _prodDecs, _pbarpTripletm1DecMap);
}



bool pbarpReaction::CheckLDrop(std::string& particleName, boost::shared_ptr<const jpcRes> theJPC){

   std::map<std::string, short> lmaxMap = pbarpEnv::instance()->lmaxParticleData();
   std::map<std::string, std::vector<short> > ldropMap = pbarpEnv::instance()->dropPbarpLForParticleData();

   auto lmaxIt = lmaxMap.find(particleName);
   if(lmaxIt != lmaxMap.end()){
      short lmaxParticle = (*lmaxIt).second;
      short lMin = 1E4;

      for(auto jpcLIt = _jpcToLMap[theJPC].begin(); jpcLIt != _jpcToLMap[theJPC].end(); ++jpcLIt){
	 if((*jpcLIt) < lMin)
	    lMin = *jpcLIt;
      }

      if(lMin > lmaxParticle)
      return false;
   }

   auto ldropIt = ldropMap.find(particleName);
   if(ldropIt != ldropMap.end()){
      bool accept = false;
      short numfound=0;
      for(auto jpcLIt = _jpcToLMap[theJPC].begin(); jpcLIt != _jpcToLMap[theJPC].end(); ++jpcLIt){
	 bool found=false;
	 for(auto ldropListIt = (*ldropIt).second.begin(); ldropListIt != (*ldropIt).second.end(); ++ldropListIt){
	    if((*ldropListIt) == (*jpcLIt)){
	       found = true;
	    }
	 }
	 if(!found)
	    accept = true;
	 else
	    numfound++;
      }

      if(accept && numfound>0){
	 Warning << "Could not drop " << numfound << " pbarp L state(s) due to JPC = " << theJPC->J << " " << theJPC->P << " " << theJPC->C <<
	    " ambiguity. Manual parameter fixing needed." << endmsg;
      }

      return accept;
   }

   return true;
}



pbarpReaction::~pbarpReaction(){
}

// void  pbarpReaction::fillMap(std::vector< boost::shared_ptr<const JPCLS> >& pbarpLSs, std::vector<boost::shared_ptr<IsobarLSDecay> >& decs, std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarLSDecay> >, pawian::Collection::SharedPtrLess > toFill){

//   std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
//   for (itJPCLS = pbarpLSs.begin(); itJPCLS != pbarpLSs.end(); ++itJPCLS){
//     std::vector<boost::shared_ptr<IsobarLSDecay> > currentIsobarVecs;

//     std::vector<boost::shared_ptr<IsobarLSDecay> >::iterator itIsobar;
//     for(itIsobar=decs.begin(); itIsobar!=decs.end(); ++itIsobar){
//       if( (*(*itIsobar)->motherJPC())==(*(*itJPCLS)) )  currentIsobarVecs.push_back(*itIsobar);
//     }
//     std::cout << "\nfill map for";
//     (*itJPCLS)->print(std::cout);
//     std::cout <<"\t with currentIsobarVecs.size()=\t" << currentIsobarVecs.size() << std::endl;
//     toFill[(*itJPCLS)]=currentIsobarVecs;   
//   }
// }

// void  pbarpReaction::fillHeliMap(std::vector< boost::shared_ptr<const JPCLS> >& pbarpLSs, std::vector<boost::shared_ptr<IsobarHeliDecay> >& decs, std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarHeliDecay> >, pawian::Collection::SharedPtrLess > toFill){

//   std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
//   for (itJPCLS = pbarpLSs.begin(); itJPCLS != pbarpLSs.end(); ++itJPCLS){
//     std::vector<boost::shared_ptr<IsobarHeliDecay> > currentHeliVecs;

//     std::vector<boost::shared_ptr<IsobarHeliDecay> >::iterator itHeli;
//     for(itHeli=decs.begin(); itHeli!=decs.end(); ++itHeli){
//       if( (*(*itHeli)->motherJPC())==(*(*itJPCLS)) )  currentHeliVecs.push_back(*itHeli);
//     }
//     std::cout << "\nfill map for";
//     (*itJPCLS)->print(std::cout);
//     std::cout <<"\t with currentHeliVecs.size()=\t" << currentHeliVecs.size() << std::endl;
//     toFill[(*itJPCLS)]=currentHeliVecs;   
//   }
// }

void pbarpReaction::print(std::ostream& os) const{
  os << "\n pbarp reaction\n";
 
  os << "\n **** JPC states for pbarp system ****** \n";
  std::vector<boost::shared_ptr<const jpcRes> >::const_iterator itJPC;
  for(itJPC=_pbarpJPCs.begin(); itJPC!=_pbarpJPCs.end(); ++itJPC){
    (*itJPC)->print(os);
  }

  os << "\n ***** decay chains *******\n";
  std::vector< boost::shared_ptr<IsobarLSDecay> >::const_iterator itIso;
  for( itIso=_prodDecs.begin(); itIso!=_prodDecs.end(); ++itIso){
    (*itIso)->print(os);
  }  
}
