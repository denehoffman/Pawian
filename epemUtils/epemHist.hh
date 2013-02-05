// epemHist class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <sstream>
#include <boost/shared_ptr.hpp>

#include "TROOT.h"


#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/AbsHist.hh"


class epemHist : public AbsHist{

public:
  epemHist(boost::shared_ptr<AbsLh>, fitParams& theFitParams);
  virtual ~epemHist();

protected:

 virtual void initRootStuff();

private:
};

