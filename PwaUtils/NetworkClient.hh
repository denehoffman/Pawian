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

#pragma once

#include <vector>

#include <boost/asio.hpp>

#include "PwaUtils/AbsChannelEnv.hh"

using boost::asio::ip::tcp;



class NetworkClient
{
  public:
   static short CLIENTMESSAGE_LOGIN;
   static short CLIENTMESSAGE_LH;
   static short CLIENTMESSAGE_HEARTBEAT;
   static short HEARTBEAT_INTERVAL; //keep this variable for timeouts

   NetworkClient(std::string serverAddress,std::string port);
   bool Login();
   bool SendLH(double llh_data, double lh_mc);
   //// heartbeat disabled
   //bool SendHeartbeat();
   bool WaitForParams();
   ChannelID channelID();
   std::vector<double>& GetParams(){return _theParams;}
   std::vector<double>& GetEventLimits(){return _eventLimits;}

  private:
   short _clientID;
   bool _paramsInitialized;
   ChannelID _channelID;
   std::vector<double> _eventLimits;
   std::string _port;
   std::string _serverAddress;
   tcp::iostream _theStream;
   //// heartbeat disabled
   //tcp::iostream _theHeartbeatStream;
   std::vector<double> _theParams;
};
