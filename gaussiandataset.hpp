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

#ifndef GAUSSIANDATASET_H
#define GAUSSIANDATASET_H

#include <Eigen/Core>
#include <string>
#include <vector>

class GaussianFCHK;

class GaussianDataset 
{
   private:
   std::string title;
   int Number_of_atoms;

   Eigen::MatrixXi atomic_numbers;
   Eigen::MatrixXd atomic_masses, 
                   cartesian_coordinates;

// hessian: force constant matrix 
//          --> hessian(i,j) = d2E/(dX_i dX_j)
//                   
// dmu:     electric dipole derivatives (APTs)
//          d mu_i / d X_k --> dmu(i,k)      
//                   
// dm:      magnetic dipole derivatives (AATs)
//          d mu_i / d Xdot_k --> dm(i,k)             
//
// NOTE: these are stored in Gaussian fchk file 
// and expressed in atomic units
//
   Eigen::MatrixXd hessian, dmu, dm;

   bool has_atomic_numbers,
        has_atomic_masses,
        has_hessian,
        has_cartesian_coordinates,
        has_electric_dipole_derivatives,
        has_magnetic_dipole_derivatives;

   public:
   GaussianDataset();
   explicit GaussianDataset(const GaussianFCHK&);

   // getters
   int GetNumberOfAtoms() const;
   Eigen::MatrixXi GetAtomicNumbers() const;
   Eigen::MatrixXd GetAtomicMasses() const;
   Eigen::MatrixXd GetCartesianCoordinates() const;
   Eigen::MatrixXd GetHessian() const;
   Eigen::MatrixXd GetElectricDipoleDerivatives() const;
   Eigen::MatrixXd GetMagneticDipoleDerivatives() const;
   
   // diagnostics
   bool Has_atomic_numbers() const;
   bool Has_atomic_masses() const;
   bool Has_hessian() const;
   bool Has_cartesian_coordinates() const;
   bool Has_electric_dipole_derivatives() const;
   bool Has_magnetic_dipole_derivatives() const;
};

#endif // GAUSSIANDATASET_H

