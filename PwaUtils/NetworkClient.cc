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

#include <iomanip>
#include <boost/asio.hpp>

#include "PwaUtils/NetworkClient.hh"
#include "PwaUtils/NetworkServer.hh"
#include "ErrLogger/ErrLogger.hh"

short NetworkClient::CLIENTMESSAGE_LOGIN = 1;
short NetworkClient::CLIENTMESSAGE_LH = 2;



NetworkClient::NetworkClient(std::string serverAddress, std::string port) :
   _port(port)
   , _serverAddress(serverAddress)
{
   _eventLimits.resize(4, 0);
   Info << "************* Client mode ****************" << endmsg;
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
   _theStream >> _eventLimits[0] >> _eventLimits[1] >> _eventLimits[2] >> _eventLimits[3];

   if(!_theStream){
      Alert << "Error: " << _theStream.error().message() << endmsg;
      return false;
   }

   Info << "Received data event range " << _eventLimits[0] << " - " << _eventLimits[1] << endmsg;
   Info << "Received mc event range " << _eventLimits[2] << " - " << _eventLimits[3] << endmsg;

   return true;
}


 
bool NetworkClient::SendLH(double llh_data, double weightSum, double lh_mc){

   _theStream.connect(_serverAddress, _port);

   if(!_theStream){
      Alert << "Could not send LH." << endmsg;
      return false;
   }

   _theStream << NetworkClient::CLIENTMESSAGE_LH << "\n";
   _theStream <<  std::setprecision(16) << llh_data << "\n" << weightSum << "\n" << lh_mc << "\n";

   return true;
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
      Alert << "Protocol error WaitForParams()" << endmsg;
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
