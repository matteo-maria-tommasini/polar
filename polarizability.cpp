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

#include "polarizability.hpp"

// support function
std::complex<double> Polarizability::alpha_ij_function(const double& dip_i, 
                                                       const double& dip_j,  
                        const double& energy, const double& omega, const double& gamma, 
                        const bool& is_electric_i, const bool& is_electric_j) const
{

   /*

      =============================================================
      NOTE ABOUT THE CHOSEN SIGN PRESCRIPTION IN THE DAMPING FACTOR
      =============================================================

      "The equal sign prescription is appropriate for the scattering situation
       when we control the initial and the final photon states. The opposite sign
       prescription is appropriate in the linear response theory when we control the
       initial state and also the form of the perturbation but we perform a summation
       over all final states. Thus only the opposite-sign convention is appropriate
       for the calculation of the atomic polarizability."

      The text above is cited from: PHYSICAL REVIEW A 76, 062106 2007 "Quantum
      electrodynamics of qubits", Iwo Bialynicki-Birula, Tomasz Sowinski, 
      DOI: 10.1103/PhysRevA.76.062106

      This agrees with:

      PHYSICAL REVIEW A, VOLUME 61, 035801 2000
      "Phenomenological damping in optical response tensors", 
      A. D. Buckingham, P. Fischer; DOI: 10.1103/PhysRevA.61.035801 

      PHYSICAL REVIEW A 76, 021803 R 2007 
      "Causal versus noncausal description of nonlinear wave mixing:
       Resolving the damping-sign controversy", 
      S. Mukamel; DOI: 10.1103/PhysRevA.76.021803 

   */

   /* 
    * ======================================
    * NOTE on transition dipoles in Gaussian
    * ======================================
    *
    * The Gaussian09 output from TDDFT calculations reports electric and
    * magnetic transition dipoles from the ground (0) to a given excited state
    * (k), i.e. mu_{0k} and m_{0k}
    *
    * Rotatory strength is defined as R_{0k} = Im(mu_{0k} m_{k0}).
    *
    * Therefore, one has to take care explicitly of the conversion from m_{0k}
    * to m_{k0}. For Hermitian operators such as the magnetic (electric)
    * dipole, this implies conjugation: m_{k0} = conj(m_{0k}).
    *
    */

   std::complex<double> c_dip_0k_i;  // <0|dipole_i|k>, where dipole could be electric/magnetic
   std::complex<double> c_dip_0k_j;  // <0|dipole_j|k>
   if (is_electric_i)
   {
      c_dip_0k_i = std::complex<double>(dip_i, 0.0);
   }
   else
   {
      c_dip_0k_i = std::complex<double>(0.0, dip_i);
   }

   if (is_electric_j)
   {
      c_dip_0k_j = std::complex<double>(dip_j, 0.0);
   }
   else
   {
      c_dip_0k_j = std::complex<double>(0.0, dip_j);
   }

   // implements the following formula:
   //
   // <0|dipole_i|k> <k|dipole_j|0> / (resonance denominator)
   //                               +
   // <k|dipole_i|0> <0|dipole_j|k> / (off-resonance denominator)
   //
   // both the electric and magnetic dipole operators are Hermitian, therefore:
   //
   // <k|dipole|0> = conj(<0|dipole|k>) 
   //
   // atomic units are assumed, i.e., energy, gamma and omega are expressed in hartree and 
   // transition dipoles in e*bohr
   //
   
   std::complex<double> c_dip_k0_i = std::conj(c_dip_0k_i);  // <k|dipole_i|0>
   std::complex<double> c_dip_k0_j = std::conj(c_dip_0k_j);  // <k|dipole_j|0>

   std::complex<double>  i_gamma_half(   0.0,  gamma / 2.0);
   std::complex<double>      c_energy(energy,          0.0);
   std::complex<double>       c_omega(omega ,          0.0);

   return(
            ( c_dip_0k_i * c_dip_k0_j ) / ( c_energy - c_omega - i_gamma_half ) 
                                        + 
            ( c_dip_k0_i * c_dip_0k_j ) / ( c_energy + c_omega + i_gamma_half ) 
         );
}

Polarizability::Polarizability(const States& states, 
                               const double& omega0, 
                               const double& gamma0) : alpha_ee(Eigen::Matrix3cd::Zero(3,3)),
                                                       alpha_em(Eigen::Matrix3cd::Zero(3,3)),
                                                       alpha_me(Eigen::Matrix3cd::Zero(3,3)),
                                                       alpha_mm(Eigen::Matrix3cd::Zero(3,3))
{

   this->omega = omega0;
   this->gamma = gamma0;

   for (unsigned int idx_state=0; idx_state<states.energy.size(); ++idx_state)
   {
      for (unsigned int i=0; i<3; ++i)
      {
         for (unsigned int j=0; j<3; ++j)
         {
            this->alpha_ee(i,j) += alpha_ij_function(states.el_dipole[idx_state](i), 
                                                     states.el_dipole[idx_state](j),
                                                     states.energy[idx_state],
                                                     this->omega, this->gamma, true, true);

            this->alpha_mm(i,j) += alpha_ij_function(states.mag_dipole[idx_state](i), 
                                                     states.mag_dipole[idx_state](j),
                                                     states.energy[idx_state],
                                                     this->omega, this->gamma, false, false);

            this->alpha_em(i,j) += alpha_ij_function(states.el_dipole[idx_state](i), 
                                                     states.mag_dipole[idx_state](j),
                                                     states.energy[idx_state],
                                                     this->omega, this->gamma, true, false);

            this->alpha_me(i,j) += alpha_ij_function(states.mag_dipole[idx_state](i), 
                                                     states.el_dipole[idx_state](j),
                                                     states.energy[idx_state],
                                                     this->omega, this->gamma, false, true);
         }
      }
   }
}

Eigen::Matrix3cd Polarizability::get_ee() const
{
   return(this->alpha_ee);
}

Eigen::Matrix3cd Polarizability::get_mm() const
{
   return(this->alpha_mm);
}

Eigen::Matrix3cd Polarizability::get_me() const
{
   return(this->alpha_me);
}

Eigen::Matrix3cd Polarizability::get_em() const
{
   return(this->alpha_em);
}

void Polarizability::Print() const
{
   Eigen::IOFormat format(4, 0, ", ", "\n", "# ", "");

   const double factor = 1.0 / 3.0;

   std::cout << "# Alpha_ee (omega = " << this->omega << "): iso = " << factor * alpha_ee.trace() << std::endl;
   std::cout << this->alpha_ee.format(format) << std::endl;

   std::cout << "# Alpha_mm (omega = " << this->omega << "): iso = " << factor * alpha_mm.trace() << std::endl;
   std::cout << this->alpha_mm.format(format) << std::endl;

   std::cout << "# Alpha_em (omega = " << this->omega << "): iso = " << factor * alpha_em.trace() << std::endl;
   std::cout << this->alpha_em.format(format) << std::endl;

   std::cout << "# Alpha_me (omega = " << this->omega << "): iso = " << factor * alpha_me.trace() << std::endl;
   std::cout << this->alpha_me.format(format) << std::endl;
}
