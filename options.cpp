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

Options::Options(const int& argc, const char *argv[])
{

  std::string option;

  if (argc == 1)
  {
     std::cout << "Run polar -h to get help." << std::endl;
     exit(0);
  }

  if (argc > 1) 
  {
     option = argv[1];
     if (option == "-h") 
     {
        std::cout << "Usage: polar [options]" << std::endl;
        std::cout << std::endl;
        std::cout << "options:      description:" << std::endl;
        std::cout << std::endl;

        std::cout << "Input file choice (controls electronic vs. vibrational polarizability calculation):" << std::endl;
        std::cout << "  -gau  file   (electronic pol.) obtain electronic states from Gaussian excited state calculation (e.g. TDDFT; default = " << this->filename_gaussian_out << ")" << std::endl;
        std::cout << std::endl;
        std::cout << "  -fchk file   (vibrational pol.) obtain vibrational states from Gaussian formatted checkpoint (default = " << this->filename_gaussian_fchk << ")" << std::endl;
        std::cout << "               NOTE: must be the checkpoint of a freq(vcd) calculation" << std::endl;
        std::cout << std::endl;
        std::cout << "  -r           carry out g-based robustness analysis (skips polarizability calculation)" << std::endl;
        std::cout << std::endl;
        std::cout << "  -d           switch on debug flag (default = " << this->debug << ")" << std::endl;
        std::cout << std::endl;

        std::cout << "Electronic polarizability:" << std::endl;
        std::cout << "  -min  E1     set minimum energy   (default = " << this->Emin    << " eV)" << std::endl;
        std::cout << "  -max  E2     set maximum energy   (default = " << this->Emax    << " eV)" << std::endl;
        std::cout << "  -g    gam    set damping to gam   (default = " << this->Gamma   << " eV)" << std::endl;
        std::cout << "  -ND    n     drop the last n states from the polarizability calculations (inactive by default)" << std::endl;
        std::cout << std::endl;

        std::cout << "Vibrational polarizability:" << std::endl;      
        std::cout << "  -wmin W1     set minimum wavenumber (default = " << this->Wmin   << " cm**-1)" << std::endl;
        std::cout << "  -wmax W2     set maximum wavenumber (default = " << this->Wmax   << " cm**-1)" << std::endl;
        std::cout << "  -wg   Wg     set damping to Wg      (default = " << this->WGamma << " cm**-1)" << std::endl;
        std::cout << std::endl;

        std::cout << "Output:" << std::endl;
        std::cout << "  -NP    N     set number of energy/wavenumber sampling points (default = " << this->Npoints << ")" << std::endl;
        std::cout << "  -dat   file  write electronic polarizability to file (default = " << this->filename_dat << ")" << std::endl;
        std::cout << "  -vdat  file  write vibrational polarizability to file (default = " << this->filename_vdat << ")" << std::endl;
        std::cout << std::endl;

        std::cout << "Compute Intensity-Carrying Modes (ICM) and save them in Molden format" << std::endl;
        std::cout << "  -icm_IR  file" << std::endl;
        std::cout << "  -icm_VCD file" << std::endl;
        std::cout << "  -icm_MAG file" << std::endl;
        std::cout << std::endl;

        exit(0);
     } 
     else 
     {

       for (unsigned int i = 0; i < argc; i++) 
       {
          option = argv[i];

          if ((option == "-min") && (argc > (i+1))) {
            this->Emin = std::atof(argv[i+1]);
            this->flag_polarizabilities = true;
          }
          if ((option == "-max") && (argc > (i+1))) {
            this->Emax = std::atof(argv[i+1]);
            this->flag_polarizabilities = true;
          }
          if ((option == "-g") && (argc > (i+1))) {
            this->Gamma = std::atof(argv[i+1]);
            this->flag_polarizabilities = true;
          }

          if ((option == "-wmin") && (argc > (i+1))) {
            this->Wmin = std::atof(argv[i+1]);
          }
          if ((option == "-wmax") && (argc > (i+1))) {
            this->Wmax = std::atof(argv[i+1]);
          }
          if ((option == "-wg") && (argc > (i+1))) {
            this->WGamma = std::atof(argv[i+1]);
          }      
          if ((option == "-NP") && (argc > (i+1))) {
            this->Npoints = std::atoi(argv[i+1]);
          }
          if ((option == "-ND") && (argc > (i+1))) {
            this->Number_of_states_to_delete = std::atoi(argv[i+1]);
            this->flag_delete_states = true;
          }
          if (option == "-d") {
            this->debug = true;
          }
          if (option == "-r") {
            this->flag_robustness = true;
          }
          if ((option == "-gau") && (argc > (i+1))) {
            this->filename_gaussian_out = argv[i+1];
            this->flag_electronic_states = true;
          }
          if ((option == "-fchk") && (argc > (i+1))) {
            this->filename_gaussian_fchk = argv[i+1];
            this->flag_vibrational_states = true;
          }
          if ((option == "-dat") && (argc > (i+1))) {
            this->filename_dat = argv[i+1];
            this->flag_dat = true;
          }
          if ((option == "-vdat") && (argc > (i+1))) {
            this->filename_vdat = argv[i+1];
            this->flag_vdat = true;
          }
          if ((option == "-icm_IR") && (argc > (i+1))) {
            this->filename_ICM = argv[i+1];
            this->flag_ICM = true;
            this->ICM_kind = "IR";
          }
          if ((option == "-icm_VCD") && (argc > (i+1))) {
            this->filename_ICM = argv[i+1];
            this->flag_ICM = true;
            this->ICM_kind = "VCD";
          }
          if ((option == "-icm_MAG") && (argc > (i+1))) {
            this->filename_ICM = argv[i+1];
            this->flag_ICM = true;
            this->ICM_kind = "MAG";
          }
       }
    }
  } 
}
