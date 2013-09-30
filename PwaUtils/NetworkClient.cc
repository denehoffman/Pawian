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

// NetworkClient class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#define _GLIBCXX_USE_NANOSLEEP

#include <iomanip>
#include <chrono>
#include <thread>
#include <math.h>
#include <limits>
#include <boost/asio.hpp>
#include <boost/bind.hpp>

#include "PwaUtils/NetworkClient.hh"
#include "PwaUtils/NetworkServer.hh"
#include "ErrLogger/ErrLogger.hh"

short NetworkClient::CLIENTMESSAGE_LOGIN = 1;
short NetworkClient::CLIENTMESSAGE_LH = 2;
short NetworkClient::CLIENTMESSAGE_HEARTBEAT = 3;
short NetworkClient::HEARTBEAT_INTERVAL = 60;

NetworkClient::NetworkClient(std::string serverAddress, std::string port) :
    _clientID(-1)
   ,_channelID(-1)
   ,_port(port)
   , _serverAddress(serverAddress)
{
   _eventLimits.resize(4, 0);
   Info << "************* Client mode ****************" << endmsg;
}



ChannelID NetworkClient::channelID(){
   if(_channelID == -1){
      Alert << "Did not receive channel id yet!" << endmsg;
      return 0;
   }
   return _channelID;
}



bool NetworkClient::Login(){

   Info << "Connecting to server " << _serverAddress << ":" << _port << endmsg;
   _theStream.connect(_serverAddress, _port);

   if(!_theStream){
      Alert << "Error: " << _theStream.error().message() << endmsg;
      return false;
   }

   _theStream << NetworkClient::CLIENTMESSAGE_LOGIN << "\n"
	      << boost::asio::ip::host_name() <<  "\n";
   _theStream >> _clientID >> _channelID
              >> _eventLimits[0] >> _eventLimits[1] >> _eventLimits[2] >> _eventLimits[3];

   if(!_theStream){
      Alert << "Error: " << _theStream.error().message() << endmsg;
      return false;
   }

   Info << "Received client id: " << _clientID << endmsg;
   Info << "Received channel id: " << _channelID << endmsg;
   Info << "Received data event range " << _eventLimits[0] << " - " << _eventLimits[1] << endmsg;
   Info << "Received mc event range " << _eventLimits[2] << " - " << _eventLimits[3] << endmsg;

   std::thread timerthread(&NetworkClient::SendHeartbeat, this);
   timerthread.detach();

   return true;
}



bool NetworkClient::SendLH(double llh_data, double lh_mc){

   _theStream.connect(_serverAddress, _port);

   if(!_theStream){
      Alert << "Could not send LH." << endmsg;
      return false;
   }

   if(isinf(llh_data) && llh_data > 0)
      llh_data=std::numeric_limits<double>::max() / 10.;
   else if(isinf(llh_data) && llh_data < 0)
      llh_data=-std::numeric_limits<double>::max() / 10.;

   _theStream << NetworkClient::CLIENTMESSAGE_LH << "\n"
              << _clientID << "\n"
              << std::setprecision(16) << llh_data << "\n" << lh_mc << "\n";

   return true;
}



bool NetworkClient::SendHeartbeat(){

   while(true){

      _theHeartbeatStream.connect(_serverAddress, _port);

      if(!_theStream){
	 Alert << "Could not send heartbeat." << endmsg;
	 return false;
      }

      _theHeartbeatStream << NetworkClient::CLIENTMESSAGE_HEARTBEAT << "\n";
      _theHeartbeatStream << _clientID << "\n";
      short answer;
      _theHeartbeatStream >> answer;

      std::this_thread::sleep_for(  std::chrono::seconds(HEARTBEAT_INTERVAL) );
   }

}



bool NetworkClient::WaitForParams(){

   if(!_theStream){
      Alert << "Stream error." << endmsg;
      return false;
   }

   short serverMessage;
   _theStream >> serverMessage;

   if(serverMessage == NetworkServer::SERVERMESSAGE_CLOSE){
      Info << "Received goodbye. Exiting." << endmsg;
      return false;
   }
   else if(serverMessage != NetworkServer::SERVERMESSAGE_PARAMS){
      Alert << "Protocol error WaitForParams() : " << serverMessage << endmsg;
      return false;
   }

   int numParams;
   _theStream >> numParams;
   _theParams.clear();

   for(int i=0; i< numParams; i++){
      double tmpParam;
      _theStream >> tmpParam;
      _theParams.push_back(tmpParam);
   }

   return true;
}
