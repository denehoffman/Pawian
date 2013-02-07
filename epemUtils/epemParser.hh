// epemParser class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include "PwaUtils/ParserBase.hh"
// Boost headers go here



class epemParser : public ParserBase 
{

  public:

  epemParser(int argc,char **argv);
  virtual ~epemParser(){;}

  const std::vector<std::string>& spinDensityNames() const { return _spinDensity;}

protected:
  virtual bool parseCommandLine(int argc,char **argv); 

  std::vector<std::string> _spinDensity;
};



