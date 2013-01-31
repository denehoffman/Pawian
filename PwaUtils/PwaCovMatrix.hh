// PwaCovMatrix class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#pragma once

#include <map>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/map.hpp>

#include "Minuit2/MnUserCovariance.h"
#include "Minuit2/MnUserParameters.h"


struct fitParams;


class PwaCovMatrix
{
  private:
   friend class boost::serialization::access;
   unsigned short _n;
   std::map<std::string, std::map<std::string, double> > _covMatrix;

  public:
   PwaCovMatrix();
   PwaCovMatrix(ROOT::Minuit2::MnUserCovariance &theMinuitCovMatrix, 
		ROOT::Minuit2::MnUserParameters &theMinuitParameters,
		fitParams &theFitParams);
   double GetElement(std::string parameter1, std::string parameter2);

   template<class Archive>
   void serialize(Archive & ar, const unsigned int version){
     ar & _n;
     ar & _covMatrix;
   }

};
