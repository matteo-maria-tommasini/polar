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

#ifndef GAUSSIANFCHK_H
#define GAUSSIANFCHK_H

#include <string>
#include <Eigen/Core>

#define NUMBER_OF_ATOMS_KEY              "Number of atoms"
#define ATOMIC_NUMBERS_KEY               "Atomic numbers" 
#define ATOMIC_MASSES_KEY                "Real atomic weights" 
#define CARTESIAN_COORDINATES_KEY        "Current cartesian coordinates"
#define HESSIAN_KEY                      "Cartesian Force Constants"
#define ELECTRIC_DIPOLE_DERIVATIVES_KEY  "Dipole Derivatives"
#define MAGNETIC_DIPOLE_DERIVATIVES_KEY  "AAT"

#define NUMBER_OF_ATOMS_MSG              "Number of atoms"
#define ATOMIC_NUMBERS_MSG               "Loaded atomic numbers" 
#define ATOMIC_MASSES_MSG                "Loaded atomic masses" 
#define CARTESIAN_COORDINATES_MSG        "Loaded Cartesian coordinates"
#define HESSIAN_MSG                      "Loaded Hessian (cartesian force constants)"
#define ELECTRIC_DIPOLE_DERIVATIVES_MSG  "Loaded electric dipole derivatives"
#define MAGNETIC_DIPOLE_DERIVATIVES_MSG  "Loaded magnetic dipole derivatives"

#define ATOMIC_NUMBERS_ERR_MSG           "NOT FOUND: atomic numbers" 
#define ATOMIC_MASSES_ERR_MSG            "NOT FOUND: atomic masses" 
#define CARTESIAN_COORDINATES_ERR_MSG    "NOT FOUND: Cartesian coordinates"
#define HESSIAN_ERR_MSG                  "NOT FOUND: Hessian (cartesian force constants)"
#define ELECTRIC_DIPOLE_DERIVATIVES_ERR_MSG "NOT FOUND: electric dipole derivatives"
#define MAGNETIC_DIPOLE_DERIVATIVES_ERR_MSG "NOT FOUND: magnetic dipole derivatives"

class GaussianFCHK
{
   private:
   std::string filename,
               title;
   int  Number_of_atoms;
   bool has_atomic_numbers,
        has_atomic_masses,
        has_hessian,
        has_cartesian_coordinates,
        has_electric_dipole_derivatives,
        has_magnetic_dipole_derivatives;

   public:
   GaussianFCHK();
   explicit GaussianFCHK(const std::string&);
   Eigen::MatrixXi ReadAtomicNumbers() const;
   Eigen::MatrixXd ReadAtomicMasses() const;
   Eigen::MatrixXd ReadCartesianCoordinates() const;
   Eigen::MatrixXd ReadHessian() const;
   Eigen::MatrixXd ReadElectricDipoleDerivatives() const;
   Eigen::MatrixXd ReadMagneticDipoleDerivatives() const;

   // getters
   int GetNumberOfAtoms() const;
   std::string GetTitle() const;

   // diagnostics
   bool Has_atomic_numbers() const;
   bool Has_atomic_masses() const;
   bool Has_hessian() const;
   bool Has_cartesian_coordinates() const;
   bool Has_electric_dipole_derivatives() const;
   bool Has_magnetic_dipole_derivatives() const;
};

#endif // GAUSSIANFCHK_H

