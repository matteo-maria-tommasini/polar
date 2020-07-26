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

#ifndef UTILS_H
#define UTILS_H

#include "states.hpp"
#include "vibstates.hpp"

void PrintStaticPolarizability(const States&);

void PrintFrequencyDependentPolarizability(const States&, 
                                           const double&, const double&, const double&, 
                                           const int&);
void WriteDatFile(const std::string&, 
                  const States&, 
                  const double&, const double&, const double&, 
                  const int&);

void WriteVibDatFile(const std::string&, 
                     const VibStates&, 
                     const double&, const double&, const double&, 
                     const int&);

#endif // UTILS_H
