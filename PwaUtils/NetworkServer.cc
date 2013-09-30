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

// Networking class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include "boost/date_time/local_time/local_time.hpp"
#include <iostream>
#include <iomanip>
#include <tuple>

#include "PwaUtils/NetworkServer.hh"
#include "PwaUtils/NetworkClient.hh"
#include "ErrLogger/ErrLogger.hh"

short NetworkServer::SERVERMESSAGE_PARAMS = 1;
short NetworkServer::SERVERMESSAGE_CLOSE = 2;
short NetworkServer::SERVERMESSAGE_OK = 3;


NetworkServer::NetworkServer(int port, unsigned short noOfClients, std::map<ChannelID, std::tuple<long, double, long> > numEventMap) :
     _port(port)
   , _clientTimeout(3*NetworkClient::HEARTBEAT_INTERVAL)
   , _globalTimeout(3*NetworkClient::HEARTBEAT_INTERVAL)
   , _noOfClients(noOfClients)
   , _closed(false)
   , _numEventMap(numEventMap)
{
   theIOService = std::shared_ptr<boost::asio::io_service>(new  boost::asio::io_service);
   theAcceptor = std::shared_ptr<tcp::acceptor>(new tcp::acceptor(*theIOService, tcp::endpoint(tcp::v4(), _port)));
   theDeadlineTimer = std::shared_ptr<boost::asio::deadline_timer>(new boost::asio::deadline_timer(*theIOService));

   for(int i=0; i<_noOfClients; i++){
      theStreams.push_back( std::shared_ptr<tcp::iostream>(new tcp::iostream) );
   }

   CalcEventDistribution(numEventMap);

   Info << "************* Server mode ****************" << endmsg;
   Info << "Listening on port " << port << endmsg;
}



bool NetworkServer::WaitForFirstClientLogin(){

   Info << "Waiting for " << _noOfClients << " clients ..." << endmsg;

   for(int i=0; i<_noOfClients; i++){

      theAcceptor->accept(*(theStreams.at(i)->rdbuf()));

      short connectionPurpose;
      *theStreams.at(i) >> connectionPurpose;

      if(connectionPurpose != NetworkClient::CLIENTMESSAGE_LOGIN){
         Alert << "ERROR: Client did not login" << endmsg;
         return false;
      }

      std::string nodeName;
      *theStreams.at(i) >> nodeName;

      Info << "Client " << nodeName << " logged in (ID " << i << ")." << endmsg;

      // Send initial client information
      *theStreams.at(i) << i << "\n"                                      // client id
                        << _eventDistribution.at(i).first << "\n"         // channel id
                        << _eventDistribution.at(i).second.at(0) << "\n"  // first data event
                        << _eventDistribution.at(i).second.at(1) << "\n"  // last data event
                        << _eventDistribution.at(i).second.at(2) << "\n"  // first mc event
                        << _eventDistribution.at(i).second.at(3) << "\n"; // last mc event

      // Store clientId->channel assignment
      _clientChannelMap[i] = _eventDistribution.at(i).first;
   }

   Info << "All clients ready." << endmsg;

   return true;
}



void NetworkServer::AcceptHandler(const boost::system::error_code& err){

   theDeadlineTimer->cancel();
   theAcceptor->cancel();
}



void NetworkServer::Timeout(const boost::system::error_code& err){

   if (err != boost::asio::error::operation_aborted){
      Alert << "Client timeout!" << endmsg;
      _closed = true;
      _globalTimeout = 9999999;
      theAcceptor->cancel();
   }
}



bool NetworkServer::WaitForLH(std::map<ChannelID, LHData>& theLHDataMap){

   for(int i=0; i<_noOfClients; i++){

      theDeadlineTimer->expires_from_now(boost::posix_time::seconds(_globalTimeout));
      theDeadlineTimer->async_wait(boost::bind(&NetworkServer::Timeout, this, boost::asio::placeholders::error ));
      theAcceptor->async_accept(*(theStreams.at(i)->rdbuf()), boost::bind(&NetworkServer::AcceptHandler, this, boost::asio::placeholders::error));
      theIOService->run();
      theIOService->reset();

      short connectionPurpose;
      *theStreams.at(i) >> connectionPurpose;

      if(connectionPurpose == NetworkClient::CLIENTMESSAGE_LOGIN){
         Alert << "Client tried to login too late." << endmsg;
         theStreams.at(i)->close();
         i--;
         continue;
      }
      else if (connectionPurpose == NetworkClient::CLIENTMESSAGE_HEARTBEAT){
         short clientID;
         *theStreams.at(i) >> clientID;

         if(!UpdateHeartbeats(clientID)){
            Timeout(boost::asio::error::timed_out);
            BroadcastClosingMessage();
         }
         else{
            *theStreams.at(i) << NetworkServer::SERVERMESSAGE_OK << "\n";
            theStreams.at(i)->flush();
            theStreams.at(i)->close();
         }

         i--;
         continue;
      }
      else if(connectionPurpose != NetworkClient::CLIENTMESSAGE_LH){
         Alert << "Protocol error in WaitForLH(): i=" << i << " CP " << connectionPurpose << endmsg;
         _closed = true;
         i--;
         return false;
      }


      double tempData, tempMc;
      short clientID;
      *theStreams.at(i) >> clientID >> tempData >> tempMc;

      ChannelID channelID = _clientChannelMap[clientID];
      theLHDataMap[channelID].logLH_data += tempData;
      theLHDataMap[channelID].LH_mc += tempMc;

      if(_closed)
      SendClosingMessage(theStreams.at(i));
   }

   return true;
}



bool NetworkServer::UpdateHeartbeats(short clientID){

   boost::posix_time::ptime now(boost::posix_time::second_clock::local_time());
   lastHeartbeats[clientID] = now;

   for(auto it = lastHeartbeats.begin(); it!= lastHeartbeats.end(); ++it){
      boost::posix_time::time_duration diff = now - (*it).second;

      if((unsigned)diff.total_seconds() >= _clientTimeout){
         Alert << "No signal from clientID " << (*it).first << " for "
               << diff.total_seconds() << " seconds." << endmsg;
         return false;
      }
   }

   return true;
}



void NetworkServer::SendParams(std::shared_ptr<tcp::iostream> destinationStream, const std::vector<double>& par){

   *destinationStream << NetworkServer::SERVERMESSAGE_PARAMS << "\n";
   *destinationStream << par.size() << "\n";

   for(auto it = par.begin(); it != par.end(); ++it){
      *destinationStream << std::setprecision(16) << *it << "\n";
   }

   destinationStream->flush();
   destinationStream->close();
}



void NetworkServer::BroadcastParams(const std::vector<double>& par){

   for(auto it = theStreams.begin(); it != theStreams.end(); ++it){
      SendParams(*it, par);
   }
}



void NetworkServer::SendClosingMessage(std::shared_ptr<tcp::iostream> destinationStream){

   *destinationStream << NetworkServer::SERVERMESSAGE_CLOSE << "\n";
   destinationStream->flush();
   destinationStream->close();
}



void NetworkServer::BroadcastClosingMessage(){

   _closed = true;
   for(auto it = theStreams.begin(); it != theStreams.end(); ++it){
      SendClosingMessage(*it);
   }
}



// This function distributes the events of the respecive channels to the
// available clients and stores the information in the _eventDistribution vector
void NetworkServer::CalcEventDistribution(std::map<ChannelID, std::tuple<long,double,long> > numEventMap){

   _eventDistribution.clear();

   // At least one client per channel is needed
   if(_noOfClients < numEventMap.size()){
      Alert << "Number of clients < number of channels!" << endmsg;
      return;
   }

   // Get total event number
   long totalEvt=0;
   for(auto it = numEventMap.begin(); it!=numEventMap.end();++it){
      long channelSum = std::get<0>((*it).second) + std::get<2>((*it).second);
      totalEvt += channelSum;
   }

   // Assign client numbers to channels
   std::vector<short> numClVec;
   for(auto it = numEventMap.begin(); it!=numEventMap.end();++it){
      long channelSum = std::get<0>((*it).second) + std::get<2>((*it).second);
      numClVec.push_back((short)(((double)channelSum / (double)totalEvt) * _noOfClients + 0.5));
   }

   // Check for #clients=0 and do correction
   for(auto it = numClVec.begin(); it!=numClVec.end();++it){
      if(*it == 0){
         // Minimum number is 1
         *it = 1;

         // Find channel with highest number of clients and decrease by one
         short max=0;
         short maxid=-1;
         int i=0;
         for(auto it2 = numClVec.begin(); it2!=numClVec.end();++it2){
            if(*it2 > max){
               max=*it2;
               maxid=i;
            }
            i++;
         }
         // Check #channel-=1 > 0
         if(numClVec.at(maxid) >= 2){
            numClVec.at(maxid)--;
         }
         else{
            Alert << "Failed to decrease client number." << endmsg;
         }
      }
   }

   // Fill event number vector
   int i=0;
   for(auto it = numEventMap.begin(); it!=numEventMap.end();++it){
      int numCl = numClVec.at(i);
      i++;

      for(int i=0;i<numCl;i++){
         // Get first and last event id for each client
         long eventStepData = (long)((double)std::get<0>((*it).second) / (double)numCl);
         long firstDataEvent = i*eventStepData;
         long lastDataEvent = (i+1)*eventStepData-1;

         long eventStepMC = (long)((double)std::get<2>((*it).second) / (double)numCl);
         long firstMCEvent = i*eventStepMC;
         long lastMCEvent = (i+1)*eventStepMC-1;

         // Last client processes last event
         if(i== (numCl - 1)){
            lastDataEvent = std::get<0>((*it).second) - 1;
            lastMCEvent   = std::get<2>((*it).second) - 1;
         }

         // Store event limits
         std::vector<long> eventRanges;
         eventRanges.push_back(firstDataEvent);
         eventRanges.push_back(lastDataEvent);
         eventRanges.push_back(firstMCEvent);
         eventRanges.push_back(lastMCEvent);
         _eventDistribution.push_back(std::pair<ChannelID, std::vector<long> >((*it).first, eventRanges));
      }
   }
}

