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

#include "vibpolarizability.hpp"
#include "units.hpp"

template <typename T> T VibPolarizability::ConvertWavenumberToVibOmegaAU(const T& wavenumber)
{
   const double conversion_factor = (2.0e0 * PI) * (SPEEDC_CGS) / (ATOMIC_FREQUENCY_UNITS);
   
   return(wavenumber * conversion_factor);
}

template <typename T> T VibPolarizability::ConvertWavenumberToVibQuantumEnergyAU(const T& wavenumber)
{
   const double conversion_factor = 4.5563e-6; 

   return(wavenumber * conversion_factor);
}

std::complex<double> VibPolarizability::ComputeElectricTransitionDipoleAU(const double& dmu_dq, 
		                                                                    const double& wavenumber)
{
   double omega_au = ConvertWavenumberToVibOmegaAU(wavenumber);

   // dmudq is e/amu**0.5 and atomic units (McWeeny) require electron mass
   // therefore the conversion 1 amu = 1 electron mass * ATMAMU
   double mu_re = std::sqrt(1.0 / (2.0 * omega_au)) * dmu_dq / std::sqrt(ATMAMU),
          mu_im = 0.0;

   return(std::complex<double>(mu_re, mu_im));
}

std::complex<double> VibPolarizability::ComputeMagneticTransitionDipoleAU(const double& dm_dqdot, 
		                                                                    const double& wavenumber)
{
   double omega_au = ConvertWavenumberToVibOmegaAU(wavenumber);

   // dmdqdot is (...)/amu**0.5 and atomic units (McWeeny) require electron mass
   // therefore the conversion 1 amu = 1 electron mass * ATMAMU
   double m_re = 0.0,
          m_im = std::sqrt(omega_au / 2.0) * dm_dqdot / std::sqrt(ATMAMU);
   
   return(std::complex<double>(m_re, m_im));
}

VibPolarizability::VibPolarizability(const VibStates& V, 
                                     const double& omega_in, 
                                     const double& gamma_in) : alpha_ee(Eigen::Matrix3cd::Zero(3,3)),
                                                               alpha_em(Eigen::Matrix3cd::Zero(3,3)),
                                                               alpha_me(Eigen::Matrix3cd::Zero(3,3)),
                                                               alpha_mm(Eigen::Matrix3cd::Zero(3,3))
{
   // init
   this->omega = omega_in;
   this->gamma = gamma_in;
  
   int Nstates                                  = V.GetNumberOfStates();
   Eigen::VectorXd              wavenumbers     = V.GetWavenumbers();
   std::vector<Eigen::Vector3d>   el_dipole_der = V.GetElectricDipoleDerivative();
   std::vector<Eigen::Vector3d>  mag_dipole_der = V.GetMagneticDipoleDerivative();
   
   std::complex<double> i_gamma_half(0.0, this->gamma / 2.0);

   // Sum Over Vibrational States
   for (unsigned int idx_state=0; idx_state<Nstates; ++idx_state)
   {
      // filter out rototranslations
      if (wavenumbers(idx_state) > 25.0)
      {

         std::complex<double> resonance_denominator = 
			      (wavenumbers(idx_state) - this->omega - i_gamma_half);

         std::complex<double> off_resonance_denominator = 
			      (wavenumbers(idx_state) + this->omega + i_gamma_half);

         std::complex<double> resonance_denominator_au =
                              ConvertWavenumberToVibQuantumEnergyAU(resonance_denominator);

         std::complex<double> off_resonance_denominator_au = 
                              ConvertWavenumberToVibQuantumEnergyAU(off_resonance_denominator);

         // (i,j) terms
         for (unsigned int i=0; i<3; ++i)
         {
            // dipole derivatives originating from (01) or (10) transitions
            std::complex<double> mu01_i = 
	              ComputeElectricTransitionDipoleAU(el_dipole_der[idx_state](i), 
					                                      wavenumbers(idx_state));
            std::complex<double> m10_i = 
	              ComputeMagneticTransitionDipoleAU(mag_dipole_der[idx_state](i), 
		            		                              wavenumbers(idx_state));
            std::complex<double> mu10_i = std::conj(mu01_i);
            std::complex<double>  m01_i = std::conj( m10_i);

            for (unsigned int j=0; j<3; ++j)
            {

               // dipole derivatives originating from (01) or (10) transitions
               std::complex<double> mu01_j = 
	                 ComputeElectricTransitionDipoleAU(el_dipole_der[idx_state](j), 
	                                                     wavenumbers(idx_state));
               std::complex<double> m10_j = 
	                 ComputeMagneticTransitionDipoleAU(mag_dipole_der[idx_state](j), 
	                                                      wavenumbers(idx_state));
               std::complex<double> mu10_j = std::conj(mu01_j);
               std::complex<double>  m01_j = std::conj( m10_j);

	       // pure electric (magnetic) polarizabilities
               this->alpha_ee(i,j) += mu01_i * mu10_j /     resonance_denominator_au;
               this->alpha_ee(i,j) += mu10_i * mu01_j / off_resonance_denominator_au;
 
               this->alpha_mm(i,j) +=  m01_i *  m10_j /     resonance_denominator_au;
               this->alpha_mm(i,j) +=  m10_i *  m01_j / off_resonance_denominator_au;

	       // mixed electric-magnetic polarizabilities 
               this->alpha_em(i,j) += mu01_i *  m10_j /     resonance_denominator_au;
               this->alpha_em(i,j) += mu10_i *  m01_j / off_resonance_denominator_au;

               this->alpha_me(i,j) +=  m01_i * mu10_j /     resonance_denominator_au;
               this->alpha_me(i,j) +=  m10_i * mu01_j / off_resonance_denominator_au;
            } 
         } 
      }

   } // main loop of SOS

}

Eigen::Matrix3cd VibPolarizability::get_ee() const
{
   return(this->alpha_ee);
}

Eigen::Matrix3cd VibPolarizability::get_mm() const
{
   return(this->alpha_mm);
}

Eigen::Matrix3cd VibPolarizability::get_me() const
{
   return(this->alpha_me);
}

Eigen::Matrix3cd VibPolarizability::get_em() const
{
   return(this->alpha_em);
}

void VibPolarizability::Print() const
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

