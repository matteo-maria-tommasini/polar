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

#include "options.hpp"
#include "utils.hpp"

#include "gaussiandataset.hpp"
#include "gaussianfchk.hpp"

#include "polarizability.hpp"
#include "vibpolarizability.hpp"
#include "dissimmetry_factors.hpp"

#include "icm.hpp"

#include "robust.hpp"
#include "vcd.hpp"

int main(const int argc, const char **argv)
{

  const Options options(argc, argv);

  ////////////////////////////////////////////////////
  // sanity check: either electronic OR vibrational //
  // polarizabilities, not both                     //
  ////////////////////////////////////////////////////
  if (options.flag_vibrational_states && options.flag_electronic_states)
  {
     std::cout << "Please choose either vibrational or electronic polarizability calculation." << std::endl;
     std::cout << "Use either the -fchk (vibrational) or -gau (electronic) option." << std::endl;
     std::cout << "Abort." << std::endl;
     exit(0);
  }

  /////////////////////////
  // robustness analysis //
  /////////////////////////
  if (options.flag_robustness)
  {
    if (options.flag_vibrational_states)
    { 
       const VibStates vibstates = VibStates(options.filename_gaussian_fchk);
       const Robust R = Robust(vibstates);
       R.Print();
       exit(0);
    }
    else if (options.flag_electronic_states) 
    {
       const States states = States(options.filename_gaussian_out);
       const Robust R = Robust(states);
       R.Print();
       exit(0);
    }
  }

  ////////////////////////////////
  // vibrational polarizability //
  ////////////////////////////////
  if (options.flag_vibrational_states)
  { 
     const VibStates vibstates = VibStates(options.filename_gaussian_fchk);
     const DissimmetryFactors g = DissimmetryFactors(vibstates);
     vibstates.PrintReport();
     g.Print();
     if (options.flag_vdat)
     {
        WriteVibDatFile(options.filename_vdat, 
                        vibstates, 
                        options.Wmin, options.Wmax, options.WGamma, options.Npoints); 
     }
     if (options.flag_ICM)
     {
        const GaussianFCHK fchk(options.filename_gaussian_fchk);
        const GaussianDataset gau = GaussianDataset(fchk);
        const ICM icm(gau, options.ICM_kind);
        // report ICMs
        icm.WriteMolden(gau, options.filename_ICM);
        icm.PrintReport();
     }
  }
  ///////////////////////////////
  // electronic polarizability //
  ///////////////////////////////
  else if (options.flag_electronic_states) 
  {
     States states = States(options.filename_gaussian_out);
     
     if (options.debug)
     {
        std::cout << "# Read " 
                  << states.GetNumberOfStates()
                  << " states." << std::endl; 
     }

     if (options.flag_delete_states)
     {
        std::cout << "# Deleting the last " 
                  << options.Number_of_states_to_delete 
                  << " states." << std::endl; 
        states.DestroyLastNstates(options.Number_of_states_to_delete);
     }

     if (options.flag_dat)
     {
        WriteDatFile(options.filename_dat, 
                     states, 
                     options.Emin, options.Emax, options.Gamma, options.Npoints); 
     }
     
     if (options.debug)
     {
        states.Print();
        std::cout << "# Kuhn-Thomas sum rule: " 
                  << states.SumRule() << std::endl;
        std::cout << "# oscillator and rotatory strengths: " << std::endl;
     }
 
     if (options.flag_polarizabilities)
     {
        if (options.debug)
        {
           // static (fully off-resonance) polarizability
           const Polarizability alpha0(states, 0.0, 0.0);
           alpha0.Print();
        }
        PrintStaticPolarizability(states);
        PrintFrequencyDependentPolarizability(states, 
           options.Emin, 
           options.Emax, 
           options.Gamma, 
           options.Npoints);
     }

  }  

  /////////
  // end //
  /////////

  return(0);
}

