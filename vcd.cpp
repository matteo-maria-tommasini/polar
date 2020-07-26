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

#include "vcd.hpp"
#include <iostream>

VCD::VCD() : intensities(), dmu_dq(), dm_dqdot(), m10_robust()
{
}

VCD::VCD(const VibraSolver& vibra, const GaussianDataset& gau)
{
   ComputeElectricDipoleDerivativesVsNormalCoordinates(vibra, gau);
   ComputeMagneticDipoleDerivativesVsNormalCoordinates(vibra, gau);

   int Nat = gau.GetNumberOfAtoms(); 
   this->intensities = Eigen::VectorXd::Zero(3*Nat);
   this->m10_robust  = Eigen::VectorXd::Zero(3*Nat);
   this->mu01_robust = Eigen::VectorXd::Zero(3*Nat);
   Eigen::MatrixXd wavenumbers = vibra.GetWavenumbers();

   // fix units and express dmu_dq, dm_dq in transition dipole moments: <0|mu|1>, <1|m|0>
   // (from atomic units to cgs units, namely esu * cm)
   for (unsigned int id_mode=0; id_mode < 3*Nat; ++id_mode)
   {
      double mu_factor = 
	                     ECHARGE * std::sqrt(NAV)
	                     * 
	                     std::sqrt(
		                              HBAR_CGS 
			                                / 
                         ( 4.0*PI*SPEEDC_CGS * std::abs(wavenumbers(id_mode)) )
	                     );
//    cgs units esu*cm
      Eigen::Vector3d mu01 = 
      { 
         this->dmu_dq[id_mode](0) * mu_factor,
         this->dmu_dq[id_mode](1) * mu_factor,
         this->dmu_dq[id_mode](2) * mu_factor
      };

//    warning: instead of Stephens' factor, we use below Buckingham's factor
//    which differs by a missing factor of 2 (and also adopts SI instead of
//    CGS).  This is because in polar code we strive to be consistent with
//    Buckingham notation, to favour the SI system and atomic units. This
//    requires multiplication by a factor of 2 of the AATs given in Gaussian09
//    fchk file -- see ReadMagneticDipoleDerivatives() in gaussianfchk class.
      double m_factor = HBAR_CGS 
                        * std::sqrt(
                                      HBAR_CGS * PI * SPEEDC_CGS 
                                      * std::abs(wavenumbers(id_mode))
                                   )
                        * ECHARGE * RBOHR_CGS 
                        * std::sqrt(NAV) / (HBAR_CGS * SPEEDC_CGS);
//    cgs units esu*cm
      Eigen::Vector3d m10 = 
      { 
         this->dm_dqdot[id_mode](0) * m_factor,
         this->dm_dqdot[id_mode](1) * m_factor,
         this->dm_dqdot[id_mode](2) * m_factor
      };
 
      // VCD intensity -- cgs units: (esu*cm)**2 
      this->intensities(id_mode) = mu01.dot(m10);
      // unit vector along the electric transition dipole
      Eigen::Vector3d u = mu01 / mu01.norm(); 
      // magnetic dipole at the robust point
      this->m10_robust(id_mode) = m10.dot(u); 
      // electric transition dipole at the robust point
      this->mu01_robust(id_mode) = mu01.norm(); 
   }
}

void VCD::ComputeElectricDipoleDerivativesVsNormalCoordinates(const VibraSolver& vibra, 
                                                              const GaussianDataset& gau)
{
   int Nat = gau.GetNumberOfAtoms(); 

   Eigen::MatrixXd  Lx = vibra.GetCartesianNuclearDisplacements();
   Eigen::MatrixXd dmu =   gau.GetElectricDipoleDerivatives();

   for (unsigned int id_mode = 0; id_mode < 3*Nat; ++id_mode)
   {
      Eigen::Vector3d dmudq = {0.0, 0.0, 0.0}; 
      for (unsigned int id_coord=0; id_coord < 3*Nat; ++id_coord)
      {
          dmudq(0) += dmu(0,id_coord) * Lx(id_coord,id_mode);
          dmudq(1) += dmu(1,id_coord) * Lx(id_coord,id_mode);
          dmudq(2) += dmu(2,id_coord) * Lx(id_coord,id_mode);
      }
      this->dmu_dq.push_back(dmudq);
   }
}

void VCD::ComputeMagneticDipoleDerivativesVsNormalCoordinates(const VibraSolver& vibra, 
                                                              const GaussianDataset& gau)
{
   int Nat = gau.GetNumberOfAtoms(); 

   Eigen::MatrixXd  Lx = vibra.GetCartesianNuclearDisplacements();
   Eigen::MatrixXd  dm =   gau.GetMagneticDipoleDerivatives();

   for (unsigned int id_mode = 0; id_mode < 3*Nat; ++id_mode)
   {
      Eigen::Vector3d mode_dm_dqdot = {0.0, 0.0, 0.0}; 
      for (unsigned int id_coord=0; id_coord < 3*Nat; ++id_coord)
      {
          mode_dm_dqdot(0) += dm(0,id_coord) * Lx(id_coord,id_mode);
          mode_dm_dqdot(1) += dm(1,id_coord) * Lx(id_coord,id_mode);
          mode_dm_dqdot(2) += dm(2,id_coord) * Lx(id_coord,id_mode);
      }
      this->dm_dqdot.push_back(mode_dm_dqdot);
   }
}

Eigen::VectorXd VCD::GetIntensities() const
{
   return(this->intensities);
}

Eigen::VectorXd VCD::Get_m10_Robust() const
{
   return(this->m10_robust);
}

Eigen::VectorXd VCD::Get_mu01_Robust() const
{
   return(this->mu01_robust);
}

std::vector<Eigen::Vector3d> VCD::GetElectricDipoleDerivativesVsNormalCoordinates() const
{
   return(this->dmu_dq);
}

std::vector<Eigen::Vector3d> VCD::GetMagneticDipoleDerivativesVsNormalCoordinates() const
{
   return(this->dm_dqdot);
}

