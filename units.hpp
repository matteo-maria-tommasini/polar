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

#define NM_PER_EV          1239.84193
#define NM_PER_WAVENUMBER  1.0e7
#define EV_PER_AU          27.21138602 

#define PI 3.1415926535897932384626433832795028841971

// physical constants in SI units (MKS, CODATA value)
#define ELECTRON_MASS       9.10938356e-31 
#define ELECTRON_CHARGE     1.6021766208e-19 
#define VACUUM_PERMITTIVITY 8.854187817e-12
#define HBAR                1.054571800e-34

// atomic units
#define HARTREE                (ELECTRON_CHARGE*ELECTRON_CHARGE) / (4.0*PI*VACUUM_PERMITTIVITY*BOHR) 
#define BOHR              4.0 * PI * VACUUM_PERMITTIVITY * HBAR * HBAR / (ELECTRON_MASS * ELECTRON_CHARGE*ELECTRON_CHARGE)
#define ALPHA_EE_UNITS    4.0 * PI * VACUUM_PERMITTIVITY * BOHR * BOHR*BOHR
#define ALPHA_EM_UNITS    4.0 * PI * VACUUM_PERMITTIVITY * HBAR * BOHR*BOHR / ELECTRON_MASS
#define E_DIPOLE_UNITS    ELECTRON_CHARGE * BOHR
#define M_DIPOLE_UNITS    ELECTRON_CHARGE * HBAR / ELECTRON_MASS

#define ATOMIC_TIME_UNITS      HBAR / (HARTREE)
#define ATOMIC_FREQUENCY_UNITS  1.0 / (ATOMIC_TIME_UNITS)

/////////////////////////////////////////////////////////////////
// these below are needed by IR and VCD intensity calculations //
// should be expressed based on the SI stuff above             //
/////////////////////////////////////////////////////////////////

// value of 1 bohr in Angstrom (10**-10 m)
#define A0_RADIUS 0.52917721092

// CGS physical constants
#define   RBOHR_CGS 1.0e-8 * A0_RADIUS 
#define ECHARGE     4.80320451e-10
#define  SPEEDC_CGS 2.99792458e10
#define    HBAR_CGS 1.054571726e-27

// Avogadro's number
#define     NAV 6.02214129e23

// atomic units constants
#define ATMTIME     2.418884326502e-17  // atomic time unit in seconds
#define ATMSPEEDC 137.035999074         // speed of light in atomic units
#define FINESTR     1.0 / 137.035999074 // fine structure constant
#define ATMAMU   1822.888479031408      // atomic mass unit expressed in electron mass


