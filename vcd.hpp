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

#ifndef VCD_H
#define VCD_H

#include "units.hpp"
#include "gaussiandataset.hpp"
#include "vibrasolver.hpp"

///////////////////////////////////////////////////////////////////////////////////
// VCD is a legacy class meant to reproduce data provided in Gaussian09 output,  //
// namely dipole and rotatory strengths (in deprecated CGS units).               //
// The simulation of VCD spectra is taken care of by the vibpolarizability class //
///////////////////////////////////////////////////////////////////////////////////
class VCD 
{
   private:
   Eigen::VectorXd intensities;
   std::vector<Eigen::Vector3d> dmu_dq, 
                                dm_dqdot;
   Eigen::VectorXd  m10_robust;
   Eigen::VectorXd mu01_robust;

   public:
   VCD();
   VCD(const VibraSolver&, const GaussianDataset&);

   // getters
   Eigen::VectorXd GetIntensities() const;
   Eigen::VectorXd Get_m10_Robust() const;
   Eigen::VectorXd Get_mu01_Robust() const;
   std::vector<Eigen::Vector3d> GetElectricDipoleDerivativesVsNormalCoordinates() const;
   std::vector<Eigen::Vector3d> GetMagneticDipoleDerivativesVsNormalCoordinates() const;

   // utility
   void ComputeElectricDipoleDerivativesVsNormalCoordinates(const VibraSolver&, const GaussianDataset&);
   void ComputeMagneticDipoleDerivativesVsNormalCoordinates(const VibraSolver&, const GaussianDataset&);

};

#endif // VCD_H

