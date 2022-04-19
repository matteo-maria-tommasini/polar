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

#include "vibstates.hpp"
#include "gaussiandataset.hpp"
#include "gaussianfchk.hpp"
#include "vibrasolver.hpp"
#include "macros.hpp"
#include <iomanip>

int VibStates::GetNumberOfStates() const
{
   return(this->wavenumbers.size());
}

std::vector<Eigen::Vector3d> VibStates::GetElectricDipoleDerivative() const
{
   return(this->dmu_dq);
}

std::vector<Eigen::Vector3d> VibStates::GetMagneticDipoleDerivative() const
{
   return(this->dm_dqdot);
}

Eigen::VectorXd VibStates::GetWavenumbers() const
{
   return(this->wavenumbers);
}

Eigen::VectorXd VibStates::Get_m10_Robust() const
{
   return(this->m10_robust);
}

Eigen::VectorXd VibStates::Get_mu01_Robust() const
{
   return(this->mu01_robust);
}

void VibStates::PrintReport() const
{
   auto N3 = this->wavenumbers.size();

   ///////////////////
   // HEADER begins //
   ///////////////////
   std::cout << BAR << std::endl;
   std::cout << "Set of " << N3 << " vibrational states (including rototranslations)" << std::endl;
   std::cout << BAR << std::endl;

   std::vector<std::string> labels;
   labels.push_back("wavenumber                                            (1/cm)");
   labels.push_back("electric dipole derivative wrt normal coordinate (x)  (e / amu**0.5)");
   labels.push_back("                                                 (y)  (e / amu**0.5)");
   labels.push_back("                                                 (z)  (e / amu**0.5)");
   labels.push_back("magnetic dipole derivative wrt normal coordinate (x)  (e**3 / (k0 * hbar amu**0.5))");
   labels.push_back("                                                 (y)  (e**3 / (k0 * hbar amu**0.5))");
   labels.push_back("                                                 (z)  (e**3 / (k0 * hbar amu**0.5))");

   auto Nlabels = labels.size();
   for (auto i=0; i<Nlabels; ++i)
   {
      std::cout << "(" << i+1 << ") " << labels[i] << std::endl; 
   }
   std::cout << BAR << std::endl;

   for (auto i=0; i<Nlabels; ++i)
   {
      std::cout << std::fixed << std::setfill(' ') 
                << std::setw(10) << std::setprecision(2)
                << i+1;
   }
   std::cout << std::endl;
   std::cout << BAR << std::endl;
   /////////////////
   // HEADER ends //
   /////////////////
  
   ///////////////////////////////////
   // print data rows, mode by mode //
   ///////////////////////////////////
   std::vector<double> values;
   for (auto i = 0; i < N3; ++i)
   {
      values.push_back(this->wavenumbers(i));

      values.push_back(this->dmu_dq[i](0));
      values.push_back(this->dmu_dq[i](1));
      values.push_back(this->dmu_dq[i](2));

      values.push_back(this->dm_dqdot[i](0));
      values.push_back(this->dm_dqdot[i](1));
      values.push_back(this->dm_dqdot[i](2));

      // final printout
      auto Nval = values.size();
      for (auto k = 0; k < Nval; ++k)
      {
         std::cout << std::fixed << std::setfill(' ') 
                   << std::setw(10) << std::setprecision(2)
                   << values[k];
      }
      std::cout << std::fixed << std::setfill(' ') << std::setw(6) 
                << "MODE" 
                << std::fixed << std::setfill(' ') << std::setw(5) 
                << (i+1)                           
                << std::endl;

      values.clear();
    }
}

VibStates::VibStates(const std::string& filename)
{
   GaussianFCHK fchk(filename);
   GaussianDataset gau = GaussianDataset(fchk);
   VibraSolver V = VibraSolver(gau);
   V.PrintReport();

   this->wavenumbers = V.GetWavenumbers();
   this->IR          = V.GetIRIntensities();
   this->VCD         = V.GetVCDIntensities();
   this->dmu_dq      = V.GetElectricDipoleDerivativesVsNormalCoordinates();
   this->dm_dqdot    = V.GetMagneticDipoleDerivativesVsNormalCoordinates();
   this->m10_robust  = V.Get_m10_Robust();
   this->mu01_robust = V.Get_mu01_Robust();
}

