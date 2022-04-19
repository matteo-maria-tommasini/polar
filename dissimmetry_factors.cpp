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

#include "dissimmetry_factors.hpp"
#include "vibpolarizability.hpp"
#include "units.hpp"
#include "macros.hpp"
#include <iomanip>

DissimmetryFactors::DissimmetryFactors(const VibStates& V)
{
   // init
   int Nstates                                  = V.GetNumberOfStates();
   this->wavenumbers                            = V.GetWavenumbers();
   std::vector<Eigen::Vector3d>   el_dipole_der = V.GetElectricDipoleDerivative(),
                                 mag_dipole_der = V.GetMagneticDipoleDerivative();

   // allocate output vectors
   this->dipole_strength = this->rotatory_strength = this->g = Eigen::VectorXd::Zero(Nstates); 
   
   // loop over Vibrational States
   for (int idx_state=0; idx_state<Nstates; ++idx_state)
   {
      Eigen::Vector3cd mu01_idx,
                       mu10_idx,
                        m01_idx,
                        m10_idx;

      // individual Cartesian components 
      for (int i=0; i<3; ++i)
      {
         // dipole derivatives originating from (01) or (10) transitions
         std::complex<double> mu01_i = 
         VibPolarizability::ComputeElectricTransitionDipoleAU(el_dipole_der[idx_state](i), 
             			                                           wavenumbers(idx_state));
         std::complex<double> m10_i = 
         VibPolarizability::ComputeMagneticTransitionDipoleAU(mag_dipole_der[idx_state](i), 
                         		                                   wavenumbers(idx_state));
         std::complex<double> mu10_i = std::conj(mu01_i);
         std::complex<double>  m01_i = std::conj( m10_i);

	      mu01_idx(i) = mu01_i;
	      mu10_idx(i) = mu10_i;
	       m10_idx(i) =  m10_i;
	       m01_idx(i) =  m01_i;
      } 
      this->dipole_strength(idx_state) = mu01_idx.dot(mu10_idx).real(); 
      this->rotatory_strength(idx_state) = mu01_idx.dot(m10_idx).imag(); 
      // this computes g in SI units, not cgs!!! 
      // that is why we have the speed of light c at the denominator (c in atomic units)
      this->g(idx_state) = 4.0*rotatory_strength(idx_state) / (ATMSPEEDC * dipole_strength(idx_state));
   } // main loop 

}

void DissimmetryFactors::Print() const
{

   ///////////////////
   // HEADER begins //
   ///////////////////
   std::cout << BAR                                                       << std::endl;
   std::cout << "DISSIMMETRY FACTORS"                                     << std::endl;
   std::cout << BAR                                                       << std::endl;

   std::vector<std::string> labels;
   labels.push_back("wavenumber                                   (1/cm)");
   labels.push_back("dipole strengths       ( D0k x 1E4,  atomic units )");
   labels.push_back("rotatory strengths     ( R0k x 1E7,  atomic units )");
   labels.push_back("dissimmetry factors    (   g x 1E4, dimensionless )");
   {
      auto N = labels.size(), i = N;
      for (i=0; i < N; ++i)
      {
         std::cout << "(" << i+1 << ") " << labels[i] << std::endl; 
      }
      std::cout << BAR << std::endl;
      for (i=0; i<N; ++i)
      {
         std::cout << std::fixed << std::setfill(' ') 
                   << std::setw(10) << std::setprecision(2)
                   << i+1;
      }
      std::cout << std::endl;
      std::cout << BAR << std::endl;
   }
   /////////////////
   // HEADER ends //
   /////////////////
   {
      auto Ng = this->g.size(), i = Ng;
      for (i = 0; i < Ng; ++i)
      {
         std::vector<double> values;

         const double d_mul = 1.0e4; 
         const double r_mul = 1.0e7; 
         const double g_mul = 1.0e4;

         values.push_back(this->wavenumbers(i));
         values.push_back(this->dipole_strength(i) * d_mul);
         values.push_back(this->rotatory_strength(i) * r_mul);
         values.push_back(this->g(i) * g_mul);

         //////////////////////////////////////////////////
         // print the data row relative to the i-th mode //
         //////////////////////////////////////////////////
         
         // print numerical data
         auto Nval = values.size(), k = Nval;
         for (k = 0; k < Nval; ++k)
         {
            std::cout << std::fixed << std::setfill(' ') 
                      << std::setw(10) << std::setprecision(2)
                      << values[k];
         }
         // end of the data row (just to be grep-friendly)
         std::cout 
                   << std::fixed << std::setfill(' ') << std::setw(6) 
                   << " (g)  MODE" 
                   << std::fixed << std::setfill(' ') << std::setw(5) 
                   << (i+1)                           
                   << std::endl;
         values.clear();
       }
   }
   

}

