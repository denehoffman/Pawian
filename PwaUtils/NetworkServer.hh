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

#include <boost/asio.hpp>

using boost::asio::ip::tcp;

class NetworkServer
{
public:
  static short SERVERMESSAGE_PARAMS;
  static short SERVERMESSAGE_CLOSE;
  
  NetworkServer(int port, short noOfClients, int numData, int numMC);
  bool WaitForLH(double& llh_data, double& weightSum, double& lh_mc);
  bool WaitForFirstClientLogin();
  void SendParams(std::shared_ptr<tcp::iostream> destinationStream, const std::vector<double>& par);
  void BroadcastParams(const std::vector<double>& par);
  void SendClosingMessage();
  int numMCs() const {return _numMC;}
  
private:

   unsigned int _port;
   unsigned int _timeout;
   unsigned short _noOfClients;	
   bool _firstLH;
   int _numData;
   int _numMC;
   std::shared_ptr<boost::asio::io_service> theIOService;
   std::shared_ptr<boost::asio::deadline_timer> theDeadlineTimer;
   std::shared_ptr<tcp::acceptor> theAcceptor;
   std::vector<std::shared_ptr<tcp::iostream>> theStreams;	 

   void Timeout(const boost::system::error_code& err);
   void AcceptHandler(const boost::system::error_code& err);
};
