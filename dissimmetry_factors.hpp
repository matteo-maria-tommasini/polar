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

#ifndef DISSIMMETRYFACTORS_H
#define DISSIMMETRYFACTORS_H

#include "vibstates.hpp"
#include <complex>
#include <cmath>

class DissimmetryFactors
{

private:
   Eigen::VectorXd  dipole_strength,
                   rotatory_strength,
                   g,
                   wavenumbers; 

public:
   explicit DissimmetryFactors(const VibStates&);
   void Print() const;

};

#endif // DISSIMMETRYFACTORS_H

