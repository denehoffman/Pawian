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
#include "FitParams/AbsParamHandler.hh"
#include "ErrLogger/ErrLogger.hh"


const short NetworkServer::SERVERMESSAGE_PARAMS = 1;
const short NetworkServer::SERVERMESSAGE_CLOSE = 2;
const short NetworkServer::SERVERMESSAGE_OK = 3;


NetworkServer::NetworkServer(int port, unsigned short noOfClients, std::map<ChannelID, std::tuple<long, double, long> >& numEventMap, std::string clientNumberWeights) :
     _port(port)
   , _clientTimeout(100*NetworkClient::HEARTBEAT_INTERVAL)
   , _globalTimeout(100*NetworkClient::HEARTBEAT_INTERVAL)
   , _noOfClients(noOfClients)
   , _noOfChannels(numEventMap.size()) 
   , _closed(false)
   , _clientParamsInitialized(false)
   , _numBroadcasted(0)
   , _clientNumberWeights(clientNumberWeights)
   , _numEventMap(numEventMap)
{
   theIOService = std::shared_ptr<boost::asio::io_service>(new  boost::asio::io_service);
   theAcceptor = std::shared_ptr<tcp::acceptor>(new tcp::acceptor(*theIOService, tcp::endpoint(tcp::v4(), _port)));
   theDeadlineTimer = std::shared_ptr<boost::asio::deadline_timer>(new boost::asio::deadline_timer(*theIOService));

   for(int i=0; i<_noOfClients; i++){
      theStreams.push_back( std::shared_ptr<tcp::iostream>(new tcp::iostream) );
   }

   InfoMsg << "************* Server mode ****************" << endmsg;
   InfoMsg << "Listening on port " << port << endmsg;

   CalcEventDistribution(numEventMap);
}



bool NetworkServer::WaitForFirstClientLogin(){

   InfoMsg << "Waiting for " << _noOfClients << " clients ..." << endmsg;

   for(int i=0; i<_noOfClients; i++){

      theAcceptor->accept(*(theStreams.at(i)->rdbuf()));

      short connectionPurpose;
      *theStreams.at(i) >> connectionPurpose;

      if(connectionPurpose == NetworkClient::CLIENTMESSAGE_HEARTBEAT){
         short clientID;
         *theStreams.at(i) >> clientID;
	 *theStreams.at(i) << NetworkServer::SERVERMESSAGE_OK << "\n";
         theStreams.at(i)->flush();
         theStreams.at(i)->close();
	 i--;
	 continue;
      }
      else if(connectionPurpose != NetworkClient::CLIENTMESSAGE_LOGIN){
         Alert << "ERROR: Client did not login. Message: " << connectionPurpose << endmsg;
	 SendClosingMessage(theStreams.at(i));
	 i--;
	 continue;
      }

      std::string nodeName;
      *theStreams.at(i) >> nodeName;

      InfoMsg << "Client " << nodeName << " logged in (ID " << i << ")." << endmsg;

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

   InfoMsg << "All clients ready." << endmsg;

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

      ChannelID channelID = _clientChannelMap.at(clientID);
      theLHDataMap[channelID].logLH_data += tempData;
      theLHDataMap[channelID].LH_mc += tempMc;

      lastLhTimes[i] = std::pair<short, boost::posix_time::ptime>(clientID, boost::posix_time::microsec_clock::local_time());

      if(_closed)
	 SendClosingMessage(theStreams.at(i));
   }

   EvalClientTiming();

   return true;
}



void NetworkServer::EvalClientTiming(){

   if(!((_numBroadcasted % 200) == 30) || lastLhTimes.size() == 0)
      return;

   boost::posix_time::ptime now(boost::posix_time::microsec_clock::local_time());

   boost::posix_time::time_duration maxdiff = now - lastLhTimes.at(0).second;

   InfoMsg << "Client timing information: " << endmsg;
   for(auto it = lastLhTimes.begin(); it!= lastLhTimes.end(); ++it){
      boost::posix_time::time_duration diff = now - (*it).second.second;

      int clientID = (*it).second.first;

      InfoMsg << "Client id " << clientID << " channel id " << _clientChannelMap.at(clientID) << " "
	   << " response time +" 
	   << std::setprecision(10) << ((double)(maxdiff.total_microseconds() - diff.total_microseconds()))/1E6 << " s" << endmsg;
   }
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



void NetworkServer::SendParams(std::shared_ptr<tcp::iostream> destinationStream, const std::vector<std::pair<unsigned int, double> >& par){

   *destinationStream << NetworkServer::SERVERMESSAGE_PARAMS << "\n";
   *destinationStream << par.size() << "\n";

   for(auto it = par.begin(); it != par.end(); ++it){
      *destinationStream << (*it).first << "\n";                           // Parameter id
      *destinationStream << std::setprecision(16) << (*it).second << "\n"; // Parameter value
   }

   destinationStream->flush();
   destinationStream->close();
}



void NetworkServer::BroadcastParams(const std::vector<double>& par){

   _numBroadcasted++;

   // The vector that will only contain changed parameters
   std::vector<std::pair<unsigned int, double> > updatedParams;

   // Send the full parameter list as initialization
   if(!_clientParamsInitialized){
      _cachedParams.resize(par.size());
      for(unsigned int i=0; i<par.size();i++){
         updatedParams.push_back(std::pair<unsigned int, double>(i, par.at(i)));
         _cachedParams.at(i) = par.at(i);
      }
      _clientParamsInitialized = true;
   }
   // Otherwise, only changed parameters
   else{
      for(unsigned int i=0; i<par.size();i++){
         if(!AbsParamHandler::CheckDoubleEquality(par.at(i), _cachedParams.at(i)) || // Check whether param has changed
            ((_numBroadcasted) % 100 == 0))                                          // But send full list every now and then
         {
            updatedParams.push_back(std::pair<unsigned int, double>(i, par.at(i)));
            _cachedParams.at(i) = par.at(i);
         }
      }
   }

   // Send changed parameters to clients
   for(auto it = theStreams.begin(); it != theStreams.end(); ++it){
      SendParams(*it, updatedParams);
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



// This function distributes the events of the respective channels to the
// available clients and stores the information in the _eventDistribution vector
void NetworkServer::CalcEventDistribution(std::map<ChannelID, std::tuple<long,double,long> >& numEventMap){

   _eventDistribution.clear();

   // At least one client per channel is needed
   if(_noOfClients < _noOfChannels){
      Alert << "Number of clients < number of channels!" << endmsg;
      return;
   }

   // Vector holding number of clients for each channel
   std::vector<short> numClVec;
   int sumCl = 0;

   // First check if number of clients (as weights) for each channel is given in the configuration file.
   // Otherwise, calc client numbers using the event numbers
   if(ReadNumClientsFromConfig(numClVec)){
      for(auto it = numClVec.begin(); it != numClVec.end(); ++it){
	 sumCl += (*it);
      }
   }
   else{
      // Get total event number
      long totalEvt=0;
      for(auto it = numEventMap.begin(); it!=numEventMap.end();++it){
	 long channelSum = std::get<0>((*it).second) + std::get<2>((*it).second);
	 totalEvt += channelSum;
      }

      // Assign client numbers to channels
      for(auto it = numEventMap.begin(); it!=numEventMap.end();++it){
	 long channelSum = std::get<0>((*it).second) + std::get<2>((*it).second);
	 short numClElem = (short)(((double)channelSum / (double)totalEvt) * _noOfClients);
	 sumCl += numClElem;
	 numClVec.push_back(numClElem);
      }
   }

   // Correct roundings
   while(sumCl != _noOfClients){
      // Find channel with lowest number of clients and increase by one
      short min=1E4;
      short minid=-1;
      int i=0;
      for(auto it = numClVec.begin(); it!=numClVec.end();++it){
	 if(*it < min){
	    min=*it;
	    minid=i;
	 }
	 i++;
      }
      
      numClVec.at(minid)++;
      sumCl++;
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

   for(unsigned int i=0; i<numClVec.size();i++){
      InfoMsg << "Number of clients for channel " << i << " : " << numClVec.at(i) << endmsg;
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

   if(_eventDistribution.size() != _noOfClients){
      Alert << "_eventDistributions has wrong size" << endmsg;
      exit(0);
   }
}



bool NetworkServer::ReadNumClientsFromConfig(std::vector<short>& numClVec){

   if(_clientNumberWeights == ""){
      return false;
   }

   std::istringstream stream(_clientNumberWeights);
   std::vector<double> weights;
   double sumOfWeights=0;

   for(int i=0; i<_noOfChannels; i++){
      double currentweight;
      if(!(stream >> currentweight)){
	 InfoMsg << "Could not read client number weights" << endmsg;
	 return false;
      }
      weights.push_back(currentweight);
      sumOfWeights += currentweight;
   }

   numClVec.clear();

   for(auto it = weights.begin(); it != weights.end(); ++it){
      numClVec.push_back((short)(_noOfClients * (*it) / sumOfWeights));
   }

   InfoMsg << "Read client number weights from configuration file." << endmsg;

   return true;
}
