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

#include "robust.hpp"
#include "vcd.hpp"
#include "vibstates.hpp"
#include "states.hpp"
#include "macros.hpp"
#include "units.hpp"
#include <iomanip>

Robust::Robust() : electric_transition_moment(), 
	               magnetic_transition_moment(), 
		           wavenumber(), 
		           transition_energy() 
{
   this->is_vibrational = false;
   this->is_electronic  = false;
}

Robust::Robust(const VibStates& S) : electric_transition_moment(S.Get_mu01_Robust()),
                                     magnetic_transition_moment(S.Get_m10_Robust()),
                                     wavenumber(S.GetWavenumbers()),
                                     transition_energy()
{
   this->is_vibrational = true;
}

Robust::Robust(const States& S) : electric_transition_moment(S.GetElectricTransitionDipoleNorm()),
                                  magnetic_transition_moment(S.GetMagneticTransitionDipoleRobustNorm()),
                                  wavenumber(),
                                  transition_energy(S.GetEnergy())                               
{
   this->is_electronic = true;
}

void Robust::Print() const
{
   unsigned int Nstates = 0;
   if (this->is_electronic)
   {
   	  Nstates = this->transition_energy.size();
   }   
   else if (this->is_vibrational)
   {
   	  Nstates = this->wavenumber.size();
   }
   
   ///////////////////
   // HEADER begins //
   ///////////////////
   std::cout << BAR << std::endl;
   if (this->is_electronic)
   {
      std::cout << "CD robustness analysis: electronic states" << std::endl;
   }
   else if (this->is_vibrational)
   {
   	  std::cout << "CD robustness analysis: vibrational states (including roto-translations)" << std::endl;
   }
   std::cout << BAR << std::endl;

   std::vector<std::string> labels;
   if (this->is_electronic)
   {
      labels.push_back("transition energy                 ( eV )");
      labels.push_back("electric transition dipole norm   ( mu, at. units )");
      labels.push_back("robust magnetic transition dipole ( m*, at. units )");
      labels.push_back("dissymmetry factor                ( g = (4 m*) / (c * mu) )");
   }
   else if (this->is_vibrational)
   {
      labels.push_back("wavenumber                        ( 1 / cm )");
      labels.push_back("electric transition dipole norm   ( mu x 1e20, in esu * cm )");
      labels.push_back("robust magnetic transition dipole ( m* x 1e24, in esu * cm )");
      labels.push_back("dissymmetry factor                (  g x 1.0E4 )");
   }
   
   for (unsigned int i=0; i<labels.size(); ++i)
   {
      std::cout << "(" << i+1 << ") " << labels[i] << std::endl; 
   }
   std::cout << BAR << std::endl;

   for (unsigned int i=0; i<labels.size(); ++i)
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
  
   /////////////////////////////////////
   // print data rows, state by state //
   /////////////////////////////////////
   std::vector<double> values;
   for (unsigned int i=0; i<Nstates; ++i)
   {  

   	  if (this->is_electronic)
      {
         values.push_back(EV_PER_AU * this->transition_energy(i));
         values.push_back(this->electric_transition_moment(i));
         values.push_back(this->magnetic_transition_moment(i));
         // SI formula for g
         values.push_back(4.0 * this->magnetic_transition_moment(i) / (ATMSPEEDC * this->electric_transition_moment(i)));
      }
      else if (this->is_vibrational)
      {
         values.push_back(this->wavenumber(i));
         values.push_back(this->electric_transition_moment(i)*1.0e20);
         values.push_back(this->magnetic_transition_moment(i)*1.0e24);
         // cgs formula for g
         values.push_back(4.0 * 1.0e4 * this->magnetic_transition_moment(i) / this->electric_transition_moment(i));
      }
      
      // final printout
      for (unsigned int k=0; k<values.size(); ++k)
      {
         std::cout << std::fixed << std::setfill(' ') 
                   << std::setw(10) << std::setprecision(3)
                   << values[k];
      }
      std::string state_label;
      if (this->is_electronic)
      {
         state_label = "STATE";
      }
      else if (this->is_vibrational)
      {
         state_label = "MODE";
      }
      std::cout << std::fixed << std::setfill(' ') << std::setw(6) 
                << state_label 
                << std::fixed << std::setfill(' ') << std::setw(5) 
                << (i+1)                           
                << std::endl;

      values.clear();
    }
}

