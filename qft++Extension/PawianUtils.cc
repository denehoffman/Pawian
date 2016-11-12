//************************************************************************//
//                                                                        //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)                      //
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

// DataUtils class definition file. -*- C++ -*-
// Copyright 2016 Bertram Kopf

#include "qft++Extension/PawianUtils.hh"

vector<LS> PawianQFT::GetValidLSWeak(const Spin &__j, const Spin &__s1, const Spin &__s2){
  vector<LS> valid_ls;
  LS ls;

  for(Spin S = abs(__s1 - __s2); S <= (__s1 + __s2); ++S){
    for(int L = (int)abs(__j - S); L <= (int)(__j + S); ++L){
        ls.L = L;
        ls.S = S;
        valid_ls.push_back(ls);
    }
  }
  return valid_ls;
}
