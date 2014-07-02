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

// NetworkServer class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#pragma once

#include <vector>
#include <memory>
#include <tuple>

#include <boost/asio.hpp>

#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"

using boost::asio::ip::tcp;

class NetworkServer
{
public:
  static const short SERVERMESSAGE_PARAMS;
  static const short SERVERMESSAGE_CLOSE;
  static const short SERVERMESSAGE_OK;

   NetworkServer(int port, unsigned short noOfClients,  std::map<ChannelID, std::tuple<long, double, long> > numEventVec, std::string clientNumberWeights);
  void CalcEventDistribution(std::map<short, std::tuple<long,double,long> > numEventMap);
  bool WaitForLH(std::map<short, LHData>& theLHDataMap);
  bool WaitForFirstClientLogin();
  void SendParams(std::shared_ptr<tcp::iostream> destinationStream, const std::vector<std::pair<unsigned int, double> >& par);
  void BroadcastParams(const std::vector<double>& par);
  void BroadcastClosingMessage();
  void SendClosingMessage(std::shared_ptr<tcp::iostream> destinationStream);
  long numMCs(ChannelID channelID) {return std::get<2>(_numEventMap[channelID]);}
  long numData(ChannelID channelID) {return std::get<0>(_numEventMap[channelID]);}
  double weightSum(ChannelID channelID) {return std::get<1>(_numEventMap[channelID]);}

private:

   unsigned int _port;
   unsigned int _clientTimeout;
   unsigned int _globalTimeout;
   unsigned short _noOfClients;
   unsigned short _noOfChannels;
   bool _closed;
   bool _clientParamsInitialized;
   int _numData;
   int _numMC;
   long _numBroadcasted;
   std::string _clientNumberWeights;
   std::shared_ptr<boost::asio::io_service> theIOService;
   std::shared_ptr<boost::asio::deadline_timer> theDeadlineTimer;
   std::shared_ptr<tcp::acceptor> theAcceptor;
   std::vector<std::shared_ptr<tcp::iostream>> theStreams;
   std::map<short, boost::posix_time::ptime > lastHeartbeats;
   std::map<short, std::pair<short, boost::posix_time::ptime > > lastLhTimes;
   std::map<short, ChannelID> _clientChannelMap;
   std::map<ChannelID, std::tuple<long, double, long> > _numEventMap;
   std::vector<std::pair<ChannelID, std::vector<long> > > _eventDistribution;
   std::vector<double> _cachedParams;

   void Timeout(const boost::system::error_code& err);
   void AcceptHandler(const boost::system::error_code& err);
   void EvalClientTiming();
   bool UpdateHeartbeats(short clientID);
   bool ReadNumClientsFromConfig(std::vector<short>& numClVec);
};
