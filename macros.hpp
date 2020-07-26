///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                 polar v1.0 - August the 15th 2019                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Copyright 2019 Matteo Tommasini                                          //
//                                                                           //
//  Licensed under the Apache License, Version 2.0 (the "License");          //
//  you may not use this file except in compliance with the License.         //
//  You may obtain a copy of the License at                                  //
//                                                                           //
//      http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                           //
//  Unless required by applicable law or agreed to in writing, software      //
//  distributed under the License is distributed on an "AS IS" BASIS,        //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
//  See the License for the specific language governing permissions and      //
//  limitations under the License.                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef MACROS_H
#define MACROS_H

// useful tensor packing-index functions (needed until Eigen will provide
// stable tensor objects) 
#define IDX2(i,j,   N)         N*i + j
#define IDX3(i,j,k, N)         N*N*i + N*j + k
#define PACKIDX(i, j)          ( i > j ? i*(i+1)/2+j : j*(j+1)/2+i )

// Levi-Civita in three dimensions 
// Note that the indexes (i,j,k) can run from 0..2 (e.g C++) or 1..3
// (e.g. Fortran) and the expression is anyway valid, since the
// expression just depends on index differences: (i-j) = ( i+1 - (j+1) ), etc.
#define EPSILON(i,j,k) (i-j)*(j-k)*(k-i)/2

// pretty printing item 
#define BAR "-------------------------------------------------------------------------------------"

#endif // MACROS_H
 

