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

#ifndef VIBSTATES_H
#define VIBSTATES_H

#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <vector>
#include "units.hpp"

class VibPolarizability;

class VibStates
{

friend class VibPolarizability;

private:
   // attributes
   Eigen::VectorXd wavenumbers,
                            IR,
                           VCD;

// electric dipole derivative wrt normal coordinates; 
// atomic units: (e / m**0.5) 
   std::vector<Eigen::Vector3d> dmu_dq;     

// magnetic dipole derivative wrt normal coordinates; 
// atomic units: (e**3 / (k0 * hbar * m**0.5))
   std::vector<Eigen::Vector3d>  dm_dqdot; 

// this comes from VCD class, where it is computed in cgs units
   Eigen::VectorXd  m10_robust;
   Eigen::VectorXd mu01_robust;

public:

   VibStates() = default;
   explicit VibStates(const std::string&);
   void PrintReport() const;

   // getters
   unsigned int GetNumberOfStates() const;
   Eigen::VectorXd              GetWavenumbers() const;
   std::vector<Eigen::Vector3d> GetElectricDipoleDerivative() const;
   std::vector<Eigen::Vector3d> GetMagneticDipoleDerivative() const;
   Eigen::VectorXd Get_m10_Robust() const;
   Eigen::VectorXd Get_mu01_Robust() const;
};

#endif // VIBSTATES_H

