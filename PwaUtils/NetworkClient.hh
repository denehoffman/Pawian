// NetworkClient class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#pragma once

#include <vector>

#include <boost/asio.hpp>

using boost::asio::ip::tcp;



class NetworkClient
{
  public:
   static short CLIENTMESSAGE_LOGIN;
   static short CLIENTMESSAGE_LH;
   NetworkClient(std::string serverAddress,std::string port);
   bool Login();
   bool SendLH(double llh_data, double weightSum, double lh_mc);
   bool WaitForParams();
   std::vector<double>& GetParams(){return _theParams;}
   std::vector<double>& GetEventLimits(){return _eventLimits;}

  private:
   std::vector<double> _eventLimits;
   std::string _port;
   std::string _serverAddress;
   tcp::iostream _theStream;
   std::vector<double> _theParams;
};
