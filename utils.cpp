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

#include <iomanip>

#include "utils.hpp"
#include "units.hpp"

#include "states.hpp"
#include "vibstates.hpp"
#include "polarizability.hpp"
#include "vibpolarizability.hpp"

void WriteVibDatFile(const std::string& filename, 
                     const VibStates& vibstates, 
                     const double& Wmin, const double& Wmax, const double& WGamma, 
                     const int& Npoints)
{

   // IO
   std::ofstream ofs;
   // setup the input file stream
   ofs.open(filename.c_str());
   if (!ofs) {
      std::cout << "Cannot write file named: " 
                << filename << std::endl;
      exit(0);
   }

   // header
   ofs << "# VIB-POLAR DAT-FILE"                                                      << std::endl;
   ofs << "#"                                                                         << std::endl;
   ofs << "# (0) GENERAL INFORMATION"                                                 << std::endl;
   ofs << "#"                                                                         << std::endl;
   ofs << "# Physical constants"                                                      << std::endl;
   ofs << "#     m                        --> electron mass"                          << std::endl;
   ofs << "#     e                        --> electron charge"                        << std::endl;
   ofs << "#     hbar                     --> Planck constant / (2*pi)"               << std::endl;
   ofs << "#     k0 = 4 * pi * epsilon_0  --> vacuum permittivity"                    << std::endl;
   ofs << "#"                                                                         << std::endl;
   ofs << "# Atomic units"                                                            << std::endl;
   ofs << "#     bohr                     --> length = k0 * hbar**2 / (m * e**2)"     << std::endl;
   ofs << "#     hartree                  --> energy = e**2 / (k0 * a0)"              << std::endl;
   ofs << "#     e * bohr                 -->  el. dipole"                            << std::endl;
   ofs << "#     e * hbar / m             --> mag. dipole"                            << std::endl;
   ofs << "#     k0 * bohr**3             -->  el.      polarizability, alpha_ee"     << std::endl;
   ofs << "#     k0 * bohr**2 * hbar / m  -->  el.-mag. polarizability, alpha_em"     << std::endl;
   ofs << "#"                                                                         << std::endl;
   ofs << "# SI values of the quantities of interest" << std::endl;
   ofs << "#     e                      = " << std::setprecision(9) << std::setw(10) << ELECTRON_CHARGE               << " " << "coulomb"                         << std::endl;
   ofs << "#     m                      = " << std::setprecision(9) << std::setw(10) << ELECTRON_MASS                 << " " << "kg"                              << std::endl;
   ofs << "#     hbar                   = " << std::setprecision(9) << std::setw(10) << HBAR                          << " " << "joules * second"                 << std::endl;
   ofs << "#     k0                     = " << std::setprecision(9) << std::setw(10) << 4.0*PI*VACUUM_PERMITTIVITY    << " " << "farad / meter"                   << std::endl;
   ofs << "#     bohr                   = " << std::setprecision(9) << std::setw(10) << BOHR                          << " " << "meters"                          << std::endl;
   ofs << "#     hartree                = " << std::setprecision(9) << std::setw(10) << HARTREE                       << " " << "joules"                          << std::endl;
   ofs << "#     atomic time            = " << std::setprecision(9) << std::setw(10) << ATOMIC_TIME_UNITS             << " " << "seconds"                         << std::endl;
   ofs << "#     atomic frequency       = " << std::setprecision(9) << std::setw(10) << ATOMIC_FREQUENCY_UNITS        << " " << "hertz"                           << std::endl;
   ofs << "#     el.  dipole      (a.u) = " << std::setprecision(9) << std::setw(10) << E_DIPOLE_UNITS                << " " << "coulomb * meters"                << std::endl;
   ofs << "#     mag. dipole      (a.u) = " << std::setprecision(9) << std::setw(10) << M_DIPOLE_UNITS                << " " << "joules / tesla"                  << std::endl;
   ofs << "#   (el.  dipole)**2   (a.u) = " << std::setprecision(9) << std::setw(10) << E_DIPOLE_UNITS*E_DIPOLE_UNITS << " " << "(coulomb * meters)**2"           << std::endl;
   ofs << "#   (el.dip)*(mag.dip) (a.u) = " << std::setprecision(9) << std::setw(10) << E_DIPOLE_UNITS*M_DIPOLE_UNITS << " " << "joule * coulomb * meter / tesla" << std::endl;
   ofs << "#     alpha_ee         (a.u) = " << std::setprecision(9) << std::setw(10) << ALPHA_EE_UNITS                << " " << "farad * meter**2"                << std::endl;
   ofs << "#     alpha_em         (a.u) = " << std::setprecision(9) << std::setw(10) << ALPHA_EM_UNITS                << " " << "coulomb * meter / tesla"         << std::endl;
   ofs << std::endl;

   // static polarizability
   {
      VibPolarizability alpha_vib_00 = VibPolarizability(vibstates, 0.0, 0.0); 
      const double factor = 1.0 / 3.0;

      // real part
      double re_ee_iso = factor * alpha_vib_00.get_ee().real().trace(), 
             re_mm_iso = factor * alpha_vib_00.get_mm().real().trace(),  
             re_em_iso = factor * alpha_vib_00.get_em().real().trace(),  
             re_me_iso = factor * alpha_vib_00.get_me().real().trace();  

      // imaginary part
      double im_ee_iso = factor * alpha_vib_00.get_ee().imag().trace(), 
             im_mm_iso = factor * alpha_vib_00.get_mm().imag().trace(),  
             im_em_iso = factor * alpha_vib_00.get_em().imag().trace(),  
             im_me_iso = factor * alpha_vib_00.get_me().imag().trace();  

      // output
      ofs << std::endl;
      ofs << "# Static polarizabilities (1/3 Trace(alpha(0,0)); atomic units - see above)"
          << std::endl;
      ofs << "# Re(ee) = " << re_ee_iso << std::endl;
      ofs << "# Im(ee) = " << im_ee_iso << std::endl;

      ofs << "# Re(em) = " << re_em_iso << std::endl;
      ofs << "# Im(em) = " << im_em_iso << std::endl;

      ofs << "# Re(me) = " << re_me_iso << std::endl;
      ofs << "# Im(me) = " << im_me_iso << std::endl;

      ofs << "# Re(mm) = " << re_mm_iso << std::endl;
      ofs << "# Im(mm) = " << im_mm_iso << std::endl;
   }

// (3) wavelength dependence
   ofs << std::endl;
   ofs << "# (3) WAVELENGTH DEPENDENCE" << std::endl;
   ofs << "#     (column 1) --> photon wavelength                 (        nm                 )" << std::endl; 
   ofs << "#     (column 2) --> photon wavenumber                 (        cm**-1             )" << std::endl; 
   ofs << "#     (column 3) --> photon frequency * 10000          ( m*e**4 / (k0**2 * hbar**3))" << std::endl; 
   ofs << "#     (column 4) --> Re(Trace(alpha_ee)) / 3           (  k0 * bohr**3             )" << std::endl; 
   ofs << "#     (column 5) --> Im(Trace(alpha_ee)) / 3           (  k0 * bohr**3             )" << std::endl; 
   ofs << "#     (column 6) --> Re(Trace(alpha_em)) / 3           (  k0 * bohr**2 * hbar / m  )" << std::endl; 
   ofs << "#     (column 7) --> Im(Trace(alpha_em)) / 3           (  k0 * bohr**2 * hbar / m  )" << std::endl; 
   ofs << "#     (column 8) --> Re(Trace(alpha_mm)) / 3           (                           )" << std::endl; 
   ofs << "#     (column 9) --> Im(Trace(alpha_mm)) / 3           (                           )" << std::endl; 

   const double DeltaW = (Wmax - Wmin) / Npoints;
   for (double W = Wmin; W < Wmax; W += DeltaW)
   {
      VibPolarizability alpha_vib = VibPolarizability(vibstates, W, WGamma); 

      // e-e part
      double alpha_vib_re_ee_iso = alpha_vib.get_ee().real().trace() / 3.0,  
             alpha_vib_im_ee_iso = alpha_vib.get_ee().imag().trace() / 3.0;  

      // e-m part
      double alpha_vib_re_em_iso = alpha_vib.get_em().real().trace() / 3.0,   
             alpha_vib_im_em_iso = alpha_vib.get_em().imag().trace() / 3.0;   

      // m-m part
      double alpha_vib_re_mm_iso = alpha_vib.get_mm().real().trace() / 3.0,   
             alpha_vib_im_mm_iso = alpha_vib.get_mm().imag().trace() / 3.0;   

      constexpr double conv = SPEEDC_CGS / (2.0 * PI * ATOMIC_FREQUENCY_UNITS);
 
      // output
      ofs 
          << std::fixed << std::setfill(' ') << std::setw( 8) << std::setprecision(2) 
          << NM_PER_WAVENUMBER / W << " " 
          << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(5) 
          << W << " " 
          << std::fixed << std::setfill(' ') << std::setw(15) << std::setprecision(5) 
          << W * conv * 10000.0 << " " 
          << std::fixed << std::setfill(' ') << std::setw(15) << std::setprecision(5) 
          << alpha_vib_re_ee_iso << " "
          << std::fixed << std::setfill(' ') << std::setw(15) << std::setprecision(5) 
          << alpha_vib_im_ee_iso << " "
          << std::fixed << std::setfill(' ') << std::setw(15) << std::setprecision(5) 
          << alpha_vib_re_em_iso << " "
          << std::fixed << std::setfill(' ') << std::setw(15) << std::setprecision(5) 
          << alpha_vib_im_em_iso << " "
          << std::fixed << std::setfill(' ') << std::setw(15) << std::setprecision(5) 
          << alpha_vib_re_mm_iso << " "
          << std::fixed << std::setfill(' ') << std::setw(15) << std::setprecision(5) 
          << alpha_vib_im_mm_iso << " +++"
          << std::endl;
   }

   // the end
   ofs.close();
}


void WriteDatFile(const std::string& filename, 
                  const States& states, 
                  const double& Emin, const double& Emax, const double& Gamma, 
                  const int& Npoints)
{

   // IO
   std::ofstream ofs;
   // setup the input file stream
   ofs.open(filename.c_str());
   if (!ofs) {
      std::cout << "Cannot write file named: " 
                << filename << std::endl;
      exit(0);
   }

   // header
   ofs << "# POLAR DAT-FILE"                                                          << std::endl;
   ofs << "#"                                                                         << std::endl;
   ofs << "# (0) GENERAL INFORMATION"                                                 << std::endl;
   ofs << "#"                                                                         << std::endl;
   ofs << "# Physical constants"                                                      << std::endl;
   ofs << "#     m                        --> electron mass"                          << std::endl;
   ofs << "#     e                        --> electron charge"                        << std::endl;
   ofs << "#     hbar                     --> Planck constant / (2*pi)"               << std::endl;
   ofs << "#     k0 = 4 * pi * epsilon_0  --> vacuum permittivity"                    << std::endl;
   ofs << "#"                                                                         << std::endl;
   ofs << "# Atomic units"                                                            << std::endl;
   ofs << "#     bohr                     --> length = k0 * hbar**2 / (m * e**2)"     << std::endl;
   ofs << "#     hartree                  --> energy = e**2 / (k0 * a0)"              << std::endl;
   ofs << "#     e * bohr                 -->  el. dipole"                            << std::endl;
   ofs << "#     e * hbar / m             --> mag. dipole"                            << std::endl;
   ofs << "#     k0 * bohr**3             -->  el.      polarizability, alpha_ee"     << std::endl;
   ofs << "#     k0 * bohr**2 * hbar / m  -->  el.-mag. polarizability, alpha_em"     << std::endl;
   ofs << "#"                                                                         << std::endl;
   ofs << "# SI values of the quantities of interest" << std::endl;
   ofs << "#     e                      = " << std::setprecision(9) << std::setw(10) << ELECTRON_CHARGE               << " " << "coulomb"                         << std::endl;
   ofs << "#     m                      = " << std::setprecision(9) << std::setw(10) << ELECTRON_MASS                 << " " << "kg"                              << std::endl;
   ofs << "#     hbar                   = " << std::setprecision(9) << std::setw(10) << HBAR                          << " " << "joules * second"                 << std::endl;
   ofs << "#     k0                     = " << std::setprecision(9) << std::setw(10) << 4.0*PI*VACUUM_PERMITTIVITY    << " " << "farad / meter"                   << std::endl;
   ofs << "#     bohr                   = " << std::setprecision(9) << std::setw(10) << BOHR                          << " " << "meters"                          << std::endl;
   ofs << "#     hartree                = " << std::setprecision(9) << std::setw(10) << HARTREE                       << " " << "joules"                          << std::endl;
   ofs << "#     el.  dipole      (a.u) = " << std::setprecision(9) << std::setw(10) << E_DIPOLE_UNITS                << " " << "coulomb * meters"                << std::endl;
   ofs << "#     mag. dipole      (a.u) = " << std::setprecision(9) << std::setw(10) << M_DIPOLE_UNITS                << " " << "joules / tesla"                  << std::endl;
   ofs << "#   (el.  dipole)**2   (a.u) = " << std::setprecision(9) << std::setw(10) << E_DIPOLE_UNITS*E_DIPOLE_UNITS << " " << "(coulomb * meters)**2"           << std::endl;
   ofs << "#   (el.dip)*(mag.dip) (a.u) = " << std::setprecision(9) << std::setw(10) << E_DIPOLE_UNITS*M_DIPOLE_UNITS << " " << "joule * coulomb * meter / tesla" << std::endl;
   ofs << "#     alpha_ee         (a.u) = " << std::setprecision(9) << std::setw(10) << ALPHA_EE_UNITS                << " " << "farad * meter**2"                << std::endl;
   ofs << "#     alpha_em         (a.u) = " << std::setprecision(9) << std::setw(10) << ALPHA_EM_UNITS                << " " << "coulomb * meter / tesla"         << std::endl;
   ofs << std::endl;

   // (1) summary of excited states properties
   ofs << "# (1) SUMMARY OF EXCITED STATE PROPERTIES" << std::endl;
   ofs << "#     (column 1)   --> state index                       (                   )" << std::endl; 
   ofs << "#     (column 2)   --> transition energy                 (       hartree     )" << std::endl; 
   ofs << "#     (column 3:5) --> transition electric dipole moment (x,y,z; e * bohr    )" << std::endl; 
   ofs << "#     (column 6:8) --> transition magnetic dipole moment (x,y,z; e * hbar / m)" << std::endl; 
   for (unsigned int idx_state = 0; idx_state < states.GetNumberOfStates(); ++idx_state)
   {
      ofs 
          << std::setfill(' ') << std::setw(4) 
          << (idx_state + 1)                          << " "
          << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(5) 
          << states.GetEnergy()[idx_state]            << " "
          << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(5) 
          << states.GetElectricDipole()[idx_state](0) << " "
          << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(5) 
          << states.GetElectricDipole()[idx_state](1) << " "
          << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(5) 
          << states.GetElectricDipole()[idx_state](2) << " "
          << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(5) 
          << states.GetMagneticDipole()[idx_state](0) << " " 
          << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(5) 
          << states.GetMagneticDipole()[idx_state](1) << " " 
          << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(5) 
          << states.GetMagneticDipole()[idx_state](2)
          << std::endl; 
   }

   // (2) oscillator picture
   ofs << std::endl;
   ofs << "# (2) OSCILLATORS PICTURE" << std::endl;
   ofs << "#     (column 1) --> state index                       (                           )" << std::endl; 
   ofs << "#     (column 2) --> transition energy                 (       hartree             )" << std::endl; 
   ofs << "#     (column 3) --> dot(mu, mu)                       ( e**2 * bohr**2            )" << std::endl; 
   ofs << "#     (column 4) --> dot( m, mu)                       ( e**2 * bohr    * hbar / m )" << std::endl; 
   for (unsigned int idx_state = 0; idx_state < states.GetNumberOfStates(); ++idx_state)
   {
      ofs 
          << std::setfill(' ') << std::setw(4) 
          << (idx_state + 1)                          << " "

          << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(5) 
          << states.GetEnergy()[idx_state]            << " "

          << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(5) 
          << states.GetElectricDipole()[idx_state].dot(states.GetElectricDipole()[idx_state]) << " "

          << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(5) 
          << states.GetElectricDipole()[idx_state].dot(states.GetMagneticDipole()[idx_state]) << " "

          << std::endl; 
   }


// (3) wavelength dependence
   ofs << std::endl;
   ofs << "# (3) WAVELENGTH DEPENDENCE" << std::endl;
   ofs << "#     (column 1) --> photon wavelength                 (        nm                 )" << std::endl; 
   ofs << "#     (column 2) --> photon energy                     (        eV                 )" << std::endl; 
   ofs << "#     (column 3) --> Re(Trace(alpha_ee)) / 3           (  k0 * bohr**3             )" << std::endl; 
   ofs << "#     (column 4) --> Im(Trace(alpha_ee)) / 3           (  k0 * bohr**3             )" << std::endl; 
   ofs << "#     (column 5) --> Re(Trace(alpha_em)) / 3           (  k0 * bohr**2 * hbar / m  )" << std::endl; 
   ofs << "#     (column 6) --> Im(Trace(alpha_em)) / 3           (  k0 * bohr**2 * hbar / m  )" << std::endl; 
   ofs << "#     (column 7) --> Re(Trace(alpha_mm)) / 3           (                           )" << std::endl; 
   ofs << "#     (column 8) --> Im(Trace(alpha_mm)) / 3           (                           )" << std::endl; 
   const double DeltaE = (Emax - Emin) / Npoints;
   for (double E = Emin; E < Emax; E += DeltaE)
   {
      Polarizability alpha = Polarizability(states, E/EV_PER_AU, Gamma/EV_PER_AU); 

      // e-e part
      double alpha_re_ee_iso = alpha.get_ee().real().trace() / 3.0,
             alpha_im_ee_iso = alpha.get_ee().imag().trace() / 3.0;  

      // e-m part
      double alpha_re_em_iso = alpha.get_em().real().trace() / 3.0,   
             alpha_im_em_iso = alpha.get_em().imag().trace() / 3.0;   

      // m-m part
      double alpha_re_mm_iso = alpha.get_mm().real().trace() / 3.0,   
             alpha_im_mm_iso = alpha.get_mm().imag().trace() / 3.0;   

      // output
      ofs 
          << std::fixed << std::setfill(' ') << std::setw( 8) << std::setprecision(2) 
          << NM_PER_EV / E << " " 
          << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(5) 
          << E << " " 
          << std::fixed << std::setfill(' ') << std::setw(15) << std::setprecision(5) 
          << alpha_re_ee_iso << " "
          << std::fixed << std::setfill(' ') << std::setw(15) << std::setprecision(5) 
          << alpha_im_ee_iso << " "
          << std::fixed << std::setfill(' ') << std::setw(15) << std::setprecision(5) 
          << alpha_re_em_iso << " "
          << std::fixed << std::setfill(' ') << std::setw(15) << std::setprecision(5) 
          << alpha_im_em_iso << " "
          << std::fixed << std::setfill(' ') << std::setw(15) << std::setprecision(5) 
          << alpha_re_mm_iso << " "
          << std::fixed << std::setfill(' ') << std::setw(15) << std::setprecision(5) 
          << alpha_im_mm_iso << " +++"
          << std::endl;
   }

   // the end
   ofs.close();
}

void PrintStaticPolarizability(const States& states)
{
   const Polarizability alpha = Polarizability(states, 0.0, 0.0); 
   const double factor = 1.0 / 3.0;

   // real part
   double alpha_re_ee_iso = factor * alpha.get_ee().real().trace(),  
          alpha_re_mm_iso = factor * alpha.get_mm().real().trace(),  
          alpha_re_em_iso = factor * alpha.get_em().real().trace(),  
          alpha_re_me_iso = factor * alpha.get_me().real().trace();  

   // imaginary part
   double alpha_im_ee_iso = factor * alpha.get_ee().imag().trace(),  
          alpha_im_mm_iso = factor * alpha.get_mm().imag().trace(),
          alpha_im_em_iso = factor * alpha.get_em().imag().trace(),  
          alpha_im_me_iso = factor * alpha.get_me().imag().trace();  

   // output
   std::cout << "# Static polarizabilities (1/3 Trace(alpha); atomic units - see above)"
             << std::endl;
   std::cout << "Re(ee) " << alpha_re_ee_iso << std::endl;
   std::cout << "Im(ee) " << alpha_im_ee_iso << std::endl;

   std::cout << "Re(em) " << alpha_re_em_iso << std::endl;
   std::cout << "Im(em) " << alpha_im_em_iso << std::endl;

   std::cout << "Re(me) " << alpha_re_me_iso << std::endl;
   std::cout << "Im(me) " << alpha_im_me_iso << std::endl;

   std::cout << "Re(mm) " << alpha_re_mm_iso << std::endl;
   std::cout << "Im(mm) " << alpha_im_mm_iso << std::endl;
}

void PrintFrequencyDependentPolarizability(const States& states, 
                                           const double& Emin, const double& Emax, const double& Gamma, 
                                           const int& Npoints)
{
   std::cout << "# iso invariants: (1/3) trace(alpha)" << std::endl; 
   std::cout << "# 1: photon energy, 2: Re(alpha_ee_iso), 3: Re(alpha_mm_iso), 4: Re(alpha_em_iso), 5: Re(alpha_me_iso)";
   std::cout << " 6: Im(alpha_ee_iso), 7: Im(alpha_mm_iso), 8: Im(alpha_em_iso), 9: Im(alpha_me_iso)" << std::endl;
   const double DeltaE = (Emax - Emin) / Npoints,
                factor = 1.0 / 3.0;
   for (double E = Emin; E < Emax; E += DeltaE)
   {
      Polarizability alpha = Polarizability(states, E/EV_PER_AU, Gamma/EV_PER_AU); 
      
      // real part
      double alpha_re_ee_iso = factor * alpha.get_ee().real().trace(),  
             alpha_re_mm_iso = factor * alpha.get_mm().real().trace(),  
             alpha_re_em_iso = factor * alpha.get_em().real().trace(),  
             alpha_re_me_iso = factor * alpha.get_me().real().trace();  
      
      // imaginary part
      double alpha_im_ee_iso = factor * alpha.get_ee().imag().trace(), 
             alpha_im_mm_iso = factor * alpha.get_mm().imag().trace(),  
             alpha_im_em_iso = factor * alpha.get_em().imag().trace(),  
             alpha_im_me_iso = factor * alpha.get_me().imag().trace();  
      
      // output
      std::cout << E << " " 
                << alpha_re_ee_iso << " "
                << alpha_re_mm_iso << " "
                << alpha_re_em_iso << " "
                << alpha_re_me_iso << " "
                << alpha_im_ee_iso << " "
                << alpha_im_mm_iso << " "
                << alpha_im_em_iso << " "
                << alpha_im_me_iso << " " 
                << std::endl;
   }
}


