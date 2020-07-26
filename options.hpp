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

#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include <iostream>

class Options
{
  public:

  double Npoints = 1200,   // sampling of the energy/wavenumber axis
         Emin    = 1.5,    // eV
         Emax    = 4.4,    // eV
         Gamma   = 0.001;  // eV

  double Wmin   =   25.0, // cm**-1
         Wmax   = 4000.0, // cm**-1
         WGamma =    5.0; // cm**-1

  bool flag_polarizabilities   = false,
       flag_electronic_states  = false,
       flag_vibrational_states = false,
       flag_ICM                = false,
       flag_robustness         = false;

  std::string filename_ICM  = "icm.mol",
                   ICM_kind = "IR";

  // controls dat file generation
  bool flag_dat = false;
  std::string filename_dat = "polar.dat";
  bool flag_vdat = false;
  std::string filename_vdat = "polar.vdat";

  // limit the number of states
  int Number_of_states_to_delete = -1;
  bool flag_delete_states = false;

  // standard input file from TDDFT calculation
  std::string filename_gaussian_out  = "gaussian.out",
              filename_gaussian_fchk = "gaussian.fchk";

  // control of extra output for debug purposes
  bool debug = false;

  Options(const int&, const char **);

};

#endif // OPTIONS_H

