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

#ifndef VIBPOLARIZABILITY_H
#define VIBPOLARIZABILITY_H

#include "vibstates.hpp"
#include <complex>
#include <cmath>

class VibPolarizability
{

private:
   Eigen::Matrix3cd alpha_ee = Eigen::Matrix3cd::Zero(3,3),
                    alpha_em = Eigen::Matrix3cd::Zero(3,3),
                    alpha_me = Eigen::Matrix3cd::Zero(3,3),
                    alpha_mm = Eigen::Matrix3cd::Zero(3,3);
   double omega = 0.0,
          gamma = 1.0;

   static std::complex<double> ComputeElectricTransitionDipoleAU(const double&, const double&);
   static std::complex<double> ComputeMagneticTransitionDipoleAU(const double&, const double&);
   template <typename T> static T ConvertWavenumberToVibOmegaAU(const T&);
   template <typename T> static T ConvertWavenumberToVibQuantumEnergyAU(const T&);

   friend class DissimmetryFactors; 

public:
   VibPolarizability(const VibStates&, const double&, const double&);
   Eigen::Matrix3cd get_ee() const;
   Eigen::Matrix3cd get_mm() const;
   Eigen::Matrix3cd get_me() const;
   Eigen::Matrix3cd get_em() const;
   void Print() const;

};

#endif // VIBPOLARIZABILITY_H

