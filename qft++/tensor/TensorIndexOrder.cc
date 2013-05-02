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
#include "qft++/tensor/TensorIndexOrder.hh"
#include <iostream>
//_____________________________________________________________________________
/** @file TensorIndexOrder.C
 *  @brief TensorIndexOrder class source file.
 */
//_____________________________________________________________________________

// void print(int _index){

//    for(int i=31; i>=0;i--){

//       int r = 1 & (_index >> i);
//       std::cout << r;

//    }
//    std::cout << std::endl;


// }


TensorIndexOrder& TensorIndexOrder::Permute() {

  if(_rank < 1) return *this;
  int level;
  bool valid = false;
  while(valid == false){
     //     print(_index);
    level = _rank - 1;
    _index += (1 << (level << 2));
    while(((*this)[level] > (int)_rank - 1) && (level >= 0)){
      this->SetIndex(level,0);
      if(level > 0){
	this->SetIndex(level-1,(*this)[level-1] + 1);
	level--; 
      }
      else{
	_index = -1;
	level--;
      }
    }
    valid = true;
    for(unsigned int i = 0; i < _rank; ++i){
      for(unsigned int j = 0; j < i; ++j){			  
	if((j < i) && (*this)[i] == (*this)[j]) valid = false;
      }
    }
  }
  return (*this); 
}
//_____________________________________________________________________________

void TensorIndexOrder::Print(ostream &__os){
  __os << "(" ;
  for(unsigned int i = 0; i < _rank; i++){
    __os << (*this)[i];
    if(i < (_rank - 1)) cout << "," ;
  }
  __os << ")" << endl;
}
//_____________________________________________________________________________
