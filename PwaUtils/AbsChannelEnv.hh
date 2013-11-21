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

// AbsChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#pragma once

#include <vector>
#include <memory>
#include "PwaUtils/AbsHist.hh"
#include "PwaUtils/WaveContribution.hh"
#include "PwaUtils/DataUtils.hh"

class Particle;
class AbsDecayList;
class IsobarDecayList;
class ParserBase;
class AbsLh;

typedef short ChannelID;


class AbsChannelEnv
{
public:

   virtual void setup(ChannelID id);
   const int noFinalStateParticles() {return _noFinalStateParticles;}
   std::vector<Particle*> finalStateParticles() {return _finalStateParticles;}
   std::vector<std::pair<Particle*, Particle*> > producedParticlePairs() {return _producedParticlePairs;}
   std::shared_ptr<AbsDecayList> absDecayList() {return _absDecList;}
   std::shared_ptr<AbsDecayList> prodDecayList() {return _prodDecList;}
   std::vector<std::vector<std::string> >& histMassSystems() {return _histMassSystems;}
   std::vector<std::shared_ptr<angleHistData> >& angleHistDataVec() {return _angleHistDataVec;}
   std::vector<std::shared_ptr<angleHistData2D> >& angleHistDataVec2D() {return _angleHistDataVec2D;}
   std::vector<std::shared_ptr<calcContributionData> >& calcContributionDataVec() {return _calcContributionDataVec;}
   bool useMassRange() {return _useMassRange;}
   bool useEvtWeight(){return _useEvtWeight;}
   double massRangeMin() {return _massMin;}
   double massRangeMax() {return _massMax;}
   std::vector<unsigned int>& particleIndicesMassRange() {return _particleIndicesMassRange;}
   Vector4<double>&  initial4Vec() {return _initial4Vec;}
   ParserBase* parser(){ return _theParser;}
   void SetLh(std::shared_ptr<AbsLh> theLh){ _theLh = theLh; }
   std::shared_ptr<AbsLh> Lh();
   ChannelID channelID(){return _channelID;}

   static short CHANNEL_PBARP;
   static short CHANNEL_EPEM;
   static short CHANNEL_RES;
   static short CHANNEL_GAMMAP;

protected:
   int _channelID;
   bool _alreadySetUp;
   bool _useEvtWeight;
   AbsChannelEnv(ParserBase* theParser);
   int _noFinalStateParticles;
   std::vector<Particle*> _finalStateParticles;
   std::vector<std::pair<Particle*, Particle*> > _producedParticlePairs;
   std::shared_ptr<AbsDecayList> _absDecList;
   std::shared_ptr<AbsDecayList> _prodDecList;
   std::vector<std::vector<std::string> > _histMassSystems;
   std::vector<std::shared_ptr<angleHistData> > _angleHistDataVec;
   std::vector<std::shared_ptr<angleHistData2D> > _angleHistDataVec2D;
   std::vector<std::shared_ptr<calcContributionData> > _calcContributionDataVec;
   std::vector<unsigned int> _particleIndicesMassRange;
   bool _useMassRange;
   double _massMin;
   double _massMax;
   Vector4<double> _initial4Vec;
   std::map<std::string, double> _preFactorMap;
   ParserBase* _theParser;
   std::shared_ptr<AbsLh> _theLh;
};
