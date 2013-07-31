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

#include "PwaUtils/NetworkServer.hh"
#include "PwaUtils/NetworkClient.hh"
#include "ErrLogger/ErrLogger.hh"

short NetworkServer::SERVERMESSAGE_PARAMS = 1;
short NetworkServer::SERVERMESSAGE_CLOSE = 2;
short NetworkServer::SERVERMESSAGE_OK = 3;


NetworkServer::NetworkServer(int port, short noOfClients, int numData, int numMC) :
     _port(port) 
   , _clientTimeout(3*NetworkClient::HEARTBEAT_INTERVAL)
   , _globalTimeout(3*NetworkClient::HEARTBEAT_INTERVAL)
   , _noOfClients(noOfClients)
   , _closed(false)
   , _numData(numData)
   , _numMC(numMC)
{
   theIOService = std::shared_ptr<boost::asio::io_service>(new  boost::asio::io_service);
   theAcceptor = std::shared_ptr<tcp::acceptor>(new tcp::acceptor(*theIOService, tcp::endpoint(tcp::v4(), _port)));
   theDeadlineTimer = std::shared_ptr<boost::asio::deadline_timer>(new boost::asio::deadline_timer(*theIOService));

   for(int i=0; i<_noOfClients; i++){
      theStreams.push_back( std::shared_ptr<tcp::iostream>(new tcp::iostream) );
   }

   Info << "************* Server mode ****************" << endmsg;
   Info << "Listening on port " << port << endmsg;
}



bool NetworkServer::WaitForFirstClientLogin(){
   
   Info << "Waiting for " << _noOfClients << " clients ..." << endmsg;

   for(int i=0; i<_noOfClients; i++){

      theAcceptor->accept(*(theStreams.at(i)->rdbuf()));

      int eventStepData = (int)((double)_numData / (double)_noOfClients);
      int firstDataEvent = i*eventStepData;
      int lastDataEvent = (i+1)*eventStepData - 1;
      
      int eventStepMC = (int)((double)_numMC / (double)_noOfClients);
      int firstMcEvent = i*eventStepMC;
      int lastMcEvent = (i+1)*eventStepMC - 1;

      if(i == (_noOfClients - 1)){
	 lastDataEvent = _numData - 1;
	 lastMcEvent = _numMC - 1;
      }

      short connectionPurpose;
      *theStreams.at(i) >> connectionPurpose;

      if(connectionPurpose != NetworkClient::CLIENTMESSAGE_LOGIN){
	 Alert << "ERROR: Client did not login" << endmsg;
	 return false;
      }

      std::string nodeName;
      *theStreams.at(i) >> nodeName;

      Info << "Client " << nodeName << " logged in." << endmsg;
      *theStreams.at(i) << firstDataEvent << "\n" << lastDataEvent << "\n" 
			<< firstMcEvent << "\n" << lastMcEvent << "\n"; 
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



bool NetworkServer::WaitForLH(double& llh_data, double& weightSum, double& lh_mc){

   for(int i=0; i<_noOfClients; i++){

      theDeadlineTimer->expires_from_now(boost::posix_time::seconds(_globalTimeout));
      theDeadlineTimer->async_wait(boost::bind(&NetworkServer::Timeout, this 
					       ,boost::asio::placeholders::error ));
      theAcceptor->async_accept(*(theStreams.at(i)->rdbuf()),
				boost::bind(&NetworkServer::AcceptHandler, this,
					    boost::asio::placeholders::error));
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
	 double clientID;
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
	 Alert << "Protocol error in WaitForLH()" << endmsg;
	 _closed = true;
	 i--;
	 return false;
      }


      double tempData, tempMc, tempWeightSum;
      *theStreams.at(i) >> tempData >> tempWeightSum >> tempMc;

      llh_data += tempData;
      lh_mc += tempMc;
      weightSum += tempWeightSum;

      if(_closed)
	 SendClosingMessage(theStreams.at(i));
   }

   return true;
}



bool NetworkServer::UpdateHeartbeats(int clientID){

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
