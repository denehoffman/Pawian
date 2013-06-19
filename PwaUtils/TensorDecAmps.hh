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

// TensorDecAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <memory>

#include "PwaUtils/AbsXdecAmp.hh"

class IsobarTensorDecay;
class AbsDecay;

class TensorDecAmps : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor 
  TensorDecAmps(std::shared_ptr<IsobarTensorDecay> theDec);
  TensorDecAmps(std::shared_ptr<AbsDecay> theDec);
  /** Destructor */
  virtual ~TensorDecAmps();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtData* theData, Spin lamFs, AbsXdecAmp* grandmaAmp);
  virtual complex<double> XdecPartAmp(Spin lamX, Spin lamDec, short fixDaughterNr,
				      EvtData* theData, Spin lamFs, AbsXdecAmp* grandmaAmp);

  virtual void print(std::ostream& os) const;
  std::vector< std::shared_ptr<const JPCLS> >& jpclsVec() {return _JPCLSs;}

  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual bool checkRecalculation(fitParams& theParamVal);
  virtual void updateFitParams(fitParams& theParamVal);

protected:
  std::vector< std::shared_ptr<const JPCLS> > _JPCLSs;
  double _factorMag;

  virtual complex<double> lsLoop(Spin lamX, EvtData* theData, Spin lam1Min, Spin lam1Max, Spin lam2Min, Spin lam2Max, bool withDecs, Spin lamFs=0 ); 
private:




};



