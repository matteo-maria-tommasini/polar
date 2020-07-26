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

#ifndef ROBUST_H
#define ROBUST_H

#include <Eigen/Core>

class VCD;
class States;
class VibStates;

class Robust 
{

private:
   // attributes
   Eigen::VectorXd electric_transition_moment, 
                   magnetic_transition_moment, 
                   wavenumber,
                   transition_energy;

   bool is_electronic  = false, 
	      is_vibrational = false;

public:
   Robust();
   explicit Robust(const VibStates&);
   explicit Robust(const States&);
   void Print() const;

};

#endif // ROBUST_H

