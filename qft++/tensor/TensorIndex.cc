// Author: Mike Williams 
/* Copyright 2008 Mike Williams (mwill@jlab.org)
 *
 * This file is part of qft++.
 *
 * qft++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * qft++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with qft++.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "qft++/tensor/TensorIndex.hh"
//_____________________________________________________________________________
/** @file TensorIndex.C
 *  @brief TensorIndex class source file.
 */
//_____________________________________________________________________________

void TensorIndex::Print(ostream &__os){
  __os << "(" ;
  for(unsigned int i = 0; i < _rank; i++){
    __os << (*this)[i];
    if(i < (_rank - 1)) cout << "," ;
  }
  __os << ")" << endl;
}
//_____________________________________________________________________________
