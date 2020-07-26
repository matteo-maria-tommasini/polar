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

#include "gaussiandataset.hpp"
#include "gaussianfchk.hpp"

// getters
int GaussianDataset::GetNumberOfAtoms() const
{
   return(this->Number_of_atoms);
}
Eigen::MatrixXd GaussianDataset::GetHessian() const
{
   return(this->hessian);
}
Eigen::MatrixXd GaussianDataset::GetAtomicMasses() const
{
   return(this->atomic_masses);
}
Eigen::MatrixXi GaussianDataset::GetAtomicNumbers() const
{
   return(this->atomic_numbers);
}
Eigen::MatrixXd GaussianDataset::GetCartesianCoordinates() const
{
   return(this->cartesian_coordinates);
}
Eigen::MatrixXd GaussianDataset::GetElectricDipoleDerivatives() const
{
   return(this->dmu);
}
Eigen::MatrixXd GaussianDataset::GetMagneticDipoleDerivatives() const
{
   return(this->dm);
}

// constructors
GaussianDataset::GaussianDataset() : title(""), Number_of_atoms(0), 
                       atomic_numbers(), 
                       atomic_masses(), 
                       cartesian_coordinates(),
                       hessian(), 
                       dmu(),
                       dm()
{
   // diagnostics
   this->has_atomic_numbers =                                                            false;
   this->has_atomic_masses =                                                             false;
   this->has_hessian =                                                                   false;
   this->has_cartesian_coordinates =                                                     false;
   this->has_electric_dipole_derivatives =                                               false;
   this->has_magnetic_dipole_derivatives =                                               false;
}

GaussianDataset::GaussianDataset(const GaussianFCHK& fchk) : title(fchk.GetTitle())
{
   this->Number_of_atoms       = fchk.GetNumberOfAtoms();
   this->atomic_numbers        = fchk.ReadAtomicNumbers();
   this->atomic_masses         = fchk.ReadAtomicMasses();
   this->cartesian_coordinates = fchk.ReadCartesianCoordinates();
   this->hessian               = fchk.ReadHessian();
   this->dmu                   = fchk.ReadElectricDipoleDerivatives();
   this->dm                    = fchk.ReadMagneticDipoleDerivatives();

   // diagnostics
   this->has_atomic_numbers = fchk.Has_atomic_numbers();
   this->has_atomic_masses = fchk.Has_atomic_masses();
   this->has_hessian = fchk.Has_hessian();
   this->has_cartesian_coordinates = fchk.Has_cartesian_coordinates();
   this->has_electric_dipole_derivatives = fchk.Has_electric_dipole_derivatives();
   this->has_magnetic_dipole_derivatives = fchk.Has_magnetic_dipole_derivatives();
}

// diagnostics
bool GaussianDataset::Has_atomic_numbers() const
{
   return(this->has_atomic_numbers);
}
bool GaussianDataset::Has_atomic_masses() const
{
   return(this->has_atomic_masses);
}
bool GaussianDataset::Has_hessian() const
{
   return(this->has_hessian);
}
bool GaussianDataset::Has_cartesian_coordinates() const
{
   return(this->has_cartesian_coordinates);
}
bool GaussianDataset::Has_electric_dipole_derivatives() const
{
   return(this->has_electric_dipole_derivatives);
}
bool GaussianDataset::Has_magnetic_dipole_derivatives() const
{
   return(this->has_magnetic_dipole_derivatives);
}

