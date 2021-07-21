//************************************************************************//
//									  //
//  Copyright 2021 Meike Kuessner (mkuessner@ep1.rub.de)                  //
//                 Bertram Kopf (bertram@ep1.rub.de)                      //  
//                 Xiaoshuai Qin (xqin@ep1.rub.de)                        //
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

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <complex>
#include <map>
#include <getopt.h>
#include <iomanip>      // std::setprecision
#include <memory>
#include <algorithm> 
#include <cmath>

#include "KMatrixExtract/TMatrixResiduePathExtr.hh"
#include "KMatrixExtract/TMatrixExtrFcn.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh" 
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceIsobar.hh"
#include "PwaDynamics/TMatrixBase.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/TMatrixNonRel.hh"
#include "PwaDynamics/KMatrixBase.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/KPoleBarrier.hh"
#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/KMatrixRelBg.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceFactory.hh"

#include "ConfigParser/KMatrixParser.hh"
#include "ConfigParser/pipiScatteringParser.hh"

#include "ErrLogger/ErrLogger.hh"
#include "Particle/PdtParser.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "FitParams/AbsPawianParamStreamer.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "FitParams/PwaCovMatrix.hh"
#include "Utils/PawianConstants.hh"

#include "ConfigParser/ParserBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/TMatrixDynamics.hh"
#include "FitParams/ParamFactory.hh"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnStrategy.h"

using namespace ROOT::Minuit2;

TMatrixResiduePathExtr::TMatrixResiduePathExtr(pipiScatteringParser* theParser) :
    TMatrixResidueExtr(theParser)
{
    _extractionMethod="Cauchy";
}

TMatrixResiduePathExtr::~TMatrixResiduePathExtr()
{
}

void TMatrixResiduePathExtr::CalcResidueAll(std::shared_ptr<AbsPawianParameters> theFitParams, 
        std::complex<double>& polePos, std::vector<ResidueProperties>& resPropCircle, 
        std::vector<ResidueProperties>& resPropEllipse, 
        std::vector<ResidueProperties>& resPropAverage){

    resPropCircle.resize(_phpVecs.size());
    resPropEllipse.resize(_phpVecs.size());
    resPropAverage.resize(_phpVecs.size());

    polePos = CalcMassWidth(theFitParams);

    InfoMsg << "\n\nm - i/2. Gamma: " << polePos.real()  << " - i/2. " << -2.*polePos.imag() << endmsg;

    for(unsigned int i=0 ; i<_phpVecs.size(); ++i) {

        ResidueProperties  currentResPropCircle;
        ResidueProperties  currentResPropEllipse;
        ResidueProperties  currentResPropAverage;

        // Circle
        std::complex<double>riemannSum_Circle=std::complex<double> (0.,0.);
        double theta = 0.;
        double dTheta = 0.001;
        double rho = 0.001;

        std::complex<double> functionArgument=std::complex<double>(0.,0.);
        std::complex<double> normFactor=std::complex<double>(0.,0.);
        std::complex<double> result_Circle =std::complex<double>(0.,0.);

        while(theta < 2.*3.141592){
            functionArgument = (polePos+ rho*std::exp(PawianConstants::i *theta));
            std::shared_ptr<TMatrixRel> currentTMat = getNewTMat();
            currentTMat->evalMatrix(functionArgument, _orbitalL);

            normFactor = std::exp(-PawianConstants::i * theta);
            riemannSum_Circle +=  1./(sqrt(_phpVecs.at(i)->factor(functionArgument, _orbitalL)) *
                    (*currentTMat)(i,i)*sqrt(_phpVecs.at(i)->factor(functionArgument, _orbitalL))) * normFactor * dTheta;
            theta += dTheta;
        }

        result_Circle = 1./(2.*3.141592*rho) * riemannSum_Circle;
        InfoMsg << "Result for circle contour: " << result_Circle << endmsg;


        // Ellipse

        std::complex<double>riemannSum_Ellipse=std::complex<double> (0.,0.);
        normFactor=std::complex<double>(0.,0.);
        theta = 0.;
        dTheta = 0.001;
        double a = 0.001;
        double b = 0.0011;

        std::complex<double> result_Ellipse =std::complex<double>(0.,0.);
        functionArgument=std::complex<double>(0.,0.);
        normFactor = std::complex<double>(0.,0.);

        while(theta < 2.*3.141592){
            functionArgument = (polePos+ std::complex<double>(a*cos(theta),b*sin(theta)));
            std::shared_ptr<TMatrixRel> currentTMat = getNewTMat();
            currentTMat->evalMatrix(functionArgument, _orbitalL);
            normFactor = 1./std::complex<double>((b*cos(theta)),(a*sin(theta)));

            riemannSum_Ellipse +=  1./(sqrt(_phpVecs.at(i)->factor(functionArgument, _orbitalL)) * 
                    (*currentTMat)(i,i)*sqrt(_phpVecs.at(i)->factor(functionArgument, _orbitalL))) * normFactor * dTheta;

            theta += dTheta;
        }

        result_Ellipse = 1./(2.*3.141592) * riemannSum_Ellipse;
        InfoMsg << "Result for elliptic contour: " << result_Ellipse << endmsg;
    
        std::complex<double> resultApprox = (result_Circle+result_Ellipse)/2.; 

        currentResPropCircle.absR=abs(1./result_Circle);
        currentResPropEllipse.absR=abs(1./result_Ellipse);
        currentResPropAverage.absR=abs(1./resultApprox);

        currentResPropCircle.theta=atan2(imag(1./result_Circle),real(1./result_Circle));
        currentResPropEllipse.theta=atan2(imag(1./result_Ellipse),real(1./result_Ellipse));
        currentResPropAverage.theta=atan2(imag(1./resultApprox),real(1./resultApprox));

        currentResPropCircle.gammai=2.*abs(1./result_Circle);
        currentResPropEllipse.gammai=2.*abs(1./result_Ellipse);
        currentResPropAverage.gammai=2.*abs(1./resultApprox);

        resPropCircle.at(i)=currentResPropCircle;
        resPropEllipse.at(i)=currentResPropEllipse;
        resPropAverage.at(i)=currentResPropAverage;

    }
}
