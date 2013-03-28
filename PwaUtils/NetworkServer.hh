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
