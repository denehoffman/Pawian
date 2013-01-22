// pbarpReaction class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "pbarpUtils/pbarpReaction.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "PwaUtils/IsobarLSDecay.hh"
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
  for(itJPC = pbarpJPCStatesAll.begin(); itJPC!=pbarpJPCStatesAll.end(); ++itJPC){
    bool acceptJPC=false;
    std::vector<std::pair<Particle*, Particle*> >::iterator itPartPairs;
    for (itPartPairs=prodPairs.begin(); itPartPairs!= prodPairs.end(); ++itPartPairs){
      std::string decName=(*itJPC)->name();
      //      std::string decName="";
      boost::shared_ptr<IsobarLSDecay> currentDec(new IsobarLSDecay( (*itJPC),itPartPairs->first, itPartPairs->second, pbarpEnv::instance(), decName));

      if (currentDec->JPCLSAmps().size()>0){
	_prodDecs.push_back(currentDec);
	acceptJPC=true;
      }
    }
    if(acceptJPC) _pbarpJPCs.push_back(*itJPC);
  }

  std::vector< boost::shared_ptr<const JPCLS> > all_JPCLSs= thepbarpStates->all_JPCLS_States();
  _pbarpJPCLSs =  extractStates(_pbarpJPCs, all_JPCLSs);

  std::vector< boost::shared_ptr<const JPCLS> > all_pbarpSingletLS = thepbarpStates->singlet_JPCLS_States();
  _pbarpJPCLSsinglet =  extractStates(_pbarpJPCs, all_pbarpSingletLS);
  fillMap(_pbarpJPCLSsinglet, _prodDecs, _pbarpSingletDecMap);

  std::vector< boost::shared_ptr<const JPCLS> > all_pbarpTriplet0LS = thepbarpStates->triplet0_JPCLS_States();
  _pbarpJPCLStriplet0 =  extractStates(_pbarpJPCs, all_pbarpTriplet0LS);
  fillMap(_pbarpJPCLStriplet0, _prodDecs, _pbarpTriplet0DecMap);

  std::vector< boost::shared_ptr<const JPCLS> > all_pbarpTripletp1LS = thepbarpStates->tripletp1_JPCLS_States();
  _pbarpJPCLStripletp1 =  extractStates(_pbarpJPCs, all_pbarpTripletp1LS);
  fillMap(_pbarpJPCLStripletp1, _prodDecs, _pbarpTripletp1DecMap);

  std::vector< boost::shared_ptr<const JPCLS> > all_pbarpTripletm1LS = thepbarpStates->tripletm1_JPCLS_States();
  _pbarpJPCLStripletm1 =  extractStates(_pbarpJPCs, all_pbarpTripletm1LS);
  fillMap(_pbarpJPCLStripletm1, _prodDecs, _pbarpTripletm1DecMap);
}

pbarpReaction::~pbarpReaction(){
}

void  pbarpReaction::fillMap(std::vector< boost::shared_ptr<const JPCLS> >& pbarpLSs, std::vector<boost::shared_ptr<IsobarLSDecay> >& decs, std::map< boost::shared_ptr<const JPCLS>, std::vector<boost::shared_ptr<IsobarLSDecay> >, pawian::Collection::SharedPtrLess > toFill){

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  for (itJPCLS = pbarpLSs.begin(); itJPCLS != pbarpLSs.end(); ++itJPCLS){
    std::vector<boost::shared_ptr<IsobarLSDecay> > currentIsobarVecs;

    std::vector<boost::shared_ptr<IsobarLSDecay> >::iterator itIsobar;
    for(itIsobar=decs.begin(); itIsobar!=decs.end(); ++itIsobar){
      if( (*(*itIsobar)->motherJPC())==(*(*itJPCLS)) )  currentIsobarVecs.push_back(*itIsobar);
    }
    std::cout << "\nfill map for";
    (*itJPCLS)->print(std::cout);
    std::cout <<"\t with currentIsobarVecs.size()=\t" << currentIsobarVecs.size() << std::endl;
    toFill[(*itJPCLS)]=currentIsobarVecs;   
  }
}


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
