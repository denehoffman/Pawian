//************************************************************************//
//                                                                        //
//  Copyright 2019 Bertram Kopf (bertram@ep1.rub.de)                      //
//                 - Ruhr-Universität Bochum                              //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

// AbsCovMatErrorExtract class definition file. -*- C++ -*-
// Copyright 2019 Bertram Kopf

#include "qaErrorExtract/AbsCovMatErrorExtract.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsLh.hh"
#include "FitParams/PwaCovMatrix.hh"
#include "ConfigParser/ParserBase.hh"
//#include "epemUtils/epemHist.hh"
//#include "pbarpUtils/pbarpHist.hh"
#include "ErrLogger/ErrLogger.hh"
//#include "FitParams/ParamFactory.hh"

#include <iostream>

AbsCovMatErrorExtract::AbsCovMatErrorExtract(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> theFitParams) :
     _calcError(false)
    , _theLh(theLh)
    , _fitParamsOriginal(theFitParams)
{
}

AbsCovMatErrorExtract::AbsCovMatErrorExtract(std::shared_ptr<AbsLh> theLh, std::shared_ptr<AbsPawianParameters> theFitParams,
				   std::shared_ptr<PwaCovMatrix> thePwaCovMatrix) :
     _calcError(true)
   , _theLh(theLh)
   , _thePwaCovMatrix(thePwaCovMatrix)
   , _fitParamsOriginal(theFitParams)
{
}

