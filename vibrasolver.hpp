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

#ifndef VIBRASOLVER_H
#define VIBRASOLVER_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues> 
#include "gaussiandataset.hpp"

class VibraSolver 
{
   private:

   Eigen::MatrixXd Lx;
   Eigen::VectorXd wavenumbers,
                   IR_intensities,
                   VCD_intensities;

   // this is computed in atomic units                
   std::vector<Eigen::Vector3d> dmu_dq,
                                 dm_dqdot;
   
   // this comes from VCD class, where it is computed in cgs units
   Eigen::VectorXd  m10_robust;
   Eigen::VectorXd mu01_robust;

   bool has_IR  = false,
        has_VCD = false;

   // utility members
   void Diagonalize(const Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::MatrixXd&);
   void SolveSecularEquation(const GaussianDataset&, Eigen::MatrixXd&, Eigen::VectorXd&);
   void WriteMatrixInOctaveFormat(const Eigen::MatrixXd&, 
                                  const std::string&, 
                                  std::ofstream&) const;

   public:

   VibraSolver() = default;
   explicit VibraSolver(const GaussianDataset&);

   // getters
   int GetNumberOfAtoms() const;
   Eigen::MatrixXd GetCartesianNuclearDisplacements() const;
   Eigen::VectorXd GetWavenumbers() const;
   Eigen::VectorXd GetIRIntensities() const;
   Eigen::VectorXd GetVCDIntensities() const;
   std::vector<Eigen::Vector3d> GetElectricDipoleDerivativesVsNormalCoordinates() const;
   std::vector<Eigen::Vector3d> GetMagneticDipoleDerivativesVsNormalCoordinates() const;
   Eigen::VectorXd Get_m10_Robust() const;
   Eigen::VectorXd Get_mu01_Robust() const;

   // output
   void PrintReport() const;

};

#endif // VIBRASOLVER_H

