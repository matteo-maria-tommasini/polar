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

#include "gaussianfchk.hpp"
#include "ioutils.hpp"
#include "macros.hpp"

std::string GaussianFCHK::GetTitle() const
{
   return(this->title);
}

int GaussianFCHK::GetNumberOfAtoms() const
{
   return(this->Number_of_atoms);
}

GaussianFCHK::GaussianFCHK() : filename(""), title("") 
{
   this->Number_of_atoms = 0;
   this->has_atomic_numbers = false;
   this->has_atomic_masses = false;
   this->has_cartesian_coordinates = false;
   this->has_hessian = false;
   this->has_electric_dipole_derivatives = false;
   this->has_magnetic_dipole_derivatives = false;
}

GaussianFCHK::GaussianFCHK(const std::string& fchk_name) : filename(fchk_name)
{
   // defaults
   this->Number_of_atoms = 0;
   this->has_atomic_numbers = false;
   this->has_atomic_masses = false;
   this->has_cartesian_coordinates = false;
   this->has_hessian = false;
   this->has_electric_dipole_derivatives = false;
   this->has_magnetic_dipole_derivatives = false;

   std::ifstream infile;
   infile.open(this->filename);

   // the first line of the formatted 
   // checkpoint file is the molecule title
   std::getline(infile, this->title);

   while (! infile.eof())
   {
      std::string line;
      std::getline(infile, line);
      if (line.find(NUMBER_OF_ATOMS_KEY) != std::string::npos) 
      {
         std::stringstream ss(extract_ints(line));
         ss >> this->Number_of_atoms;
         std::cout << NUMBER_OF_ATOMS_MSG 
                   << ": " << this->Number_of_atoms << std::endl;
      }
      if (line.find(ATOMIC_NUMBERS_KEY) != std::string::npos) 
      {
         this->has_atomic_numbers = true;
      }
      if (line.find(ATOMIC_MASSES_KEY) != std::string::npos) 
      {
         this->has_atomic_masses = true;
      }
      if (line.find(CARTESIAN_COORDINATES_KEY) != std::string::npos) 
      {
         this->has_cartesian_coordinates = true;
      }
      if (line.find(HESSIAN_KEY) != std::string::npos) 
      {
         this->has_hessian = true;
      }
      if (line.find(ELECTRIC_DIPOLE_DERIVATIVES_KEY) != std::string::npos) 
      {
         this->has_electric_dipole_derivatives = true;
      }
      if (line.find(MAGNETIC_DIPOLE_DERIVATIVES_KEY) != std::string::npos) 
      {
         this->has_magnetic_dipole_derivatives = true;
      }
   }
   infile.close();
}



Eigen::MatrixXd GaussianFCHK::ReadHessian() const
{
   if (this->has_hessian)
   {
      std::ifstream infile;
      infile.open(this->filename);
      Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3*this->Number_of_atoms,
                                                3*this->Number_of_atoms);

      std::string line;
      while (! infile.eof())
      {
         std::getline(infile, line);
         if (line.find(HESSIAN_KEY) != std::string::npos) 
         {
            int number_of_items = (3*this->Number_of_atoms)*(3*this->Number_of_atoms + 1)/2;

            std::vector<double> hessian_v = readitems<double>(infile, number_of_items);

            for (unsigned int i=0; i<3*this->Number_of_atoms; ++i) 
            {
               for (unsigned int j=i; j<3*this->Number_of_atoms; ++j) 
               {
                  H(i,j) = hessian_v[PACKIDX(i,j)];
                  H(j,i) = H(i,j); 
               }
            }
            std::cout << HESSIAN_MSG << std::endl;
         }
      }
      infile.close();
      return(H);
   }
   else
   {
      std::cout << HESSIAN_ERR_MSG << std::endl;
      Eigen::MatrixXd H; // this is empty
      return(H);
   }
}


Eigen::MatrixXi GaussianFCHK::ReadAtomicNumbers() const
{
   if (this->has_atomic_numbers)
   {
      std::ifstream infile;
      infile.open(this->filename);
      const int N = this->Number_of_atoms;
      Eigen::MatrixXi Z = Eigen::MatrixXi::Zero(N,1);

      std::string line;
      while (! infile.eof())
      {
         std::getline(infile, line);
         if (line.find(ATOMIC_NUMBERS_KEY) != std::string::npos) 
         {
            std::vector<double> Z_v = readitems<double>(infile, N);
            for (unsigned int i=0; i<N; ++i) Z(i) = Z_v[i];
            std::cout << ATOMIC_NUMBERS_MSG << std::endl;
         }
      }
      infile.close();
      return(Z);
   }
   else
   {
      std::cout << ATOMIC_NUMBERS_ERR_MSG << std::endl;
      Eigen::MatrixXi Z; // this is empty
      return(Z);
   }
}

Eigen::MatrixXd GaussianFCHK::ReadAtomicMasses() const
{
   if (this->has_atomic_masses)
   {
      std::ifstream infile;
      infile.open(this->filename);
      const int N = this->Number_of_atoms;
      Eigen::MatrixXd masses = Eigen::MatrixXd::Zero(N,1);

      std::string line;
      while (! infile.eof())
      {
         std::getline(infile, line);
         if (line.find(ATOMIC_MASSES_KEY) != std::string::npos) 
         {
            std::vector<double> masses_v = readitems<double>(infile, N);
            for (unsigned int i=0; i<N; ++i) masses(i) = masses_v[i];
            std::cout << ATOMIC_MASSES_MSG << std::endl;
         }
      }
      infile.close();
      return(masses);
   }
   else
   {
      std::cout << ATOMIC_MASSES_ERR_MSG << std::endl;
      Eigen::MatrixXd masses; // this is empty
      return(masses);
   }
}


Eigen::MatrixXd GaussianFCHK::ReadCartesianCoordinates() const
{
   if (this->has_cartesian_coordinates)
   {
      std::ifstream infile;
      infile.open(this->filename);
      const int N = this->Number_of_atoms;
      Eigen::MatrixXd X = Eigen::MatrixXd::Zero(3,N);

      std::string line;
      while (! infile.eof())
      {
         std::getline(infile, line);
         if (line.find(CARTESIAN_COORDINATES_KEY) != std::string::npos) 
         {
            std::vector<double> coord_v = readitems<double>(infile, 3*N);
            for (unsigned int i=0; i < N; ++i) 
            {
               X(0,i) = coord_v[3*i + 0];
               X(1,i) = coord_v[3*i + 1];
               X(2,i) = coord_v[3*i + 2];
            }
            std::cout << CARTESIAN_COORDINATES_MSG << std::endl;
         }
      }
      infile.close();
      return(X);
   }
   else
   {
      std::cout << CARTESIAN_COORDINATES_ERR_MSG << std::endl;
      Eigen::MatrixXd X;
      return(X);
   }
}

Eigen::MatrixXd GaussianFCHK::ReadElectricDipoleDerivatives() const
{
   if (this->has_electric_dipole_derivatives)
   {
      std::ifstream infile;
      infile.open(this->filename);
      const int N = this->Number_of_atoms;
      Eigen::MatrixXd dmu = Eigen::MatrixXd::Zero(3,3*N);

      std::string line;
      while (! infile.eof())
      {
         std::getline(infile, line);
         if (line.find(ELECTRIC_DIPOLE_DERIVATIVES_KEY) != std::string::npos) 
         {
            std::vector<double> dmu_v = readitems<double>(infile, 3*(3*N));
            for (unsigned int i=0; i < 3*N; ++i) 
            {
               dmu(0,i) = dmu_v[3*i + 0];
               dmu(1,i) = dmu_v[3*i + 1];
               dmu(2,i) = dmu_v[3*i + 2];
            }
            std::cout << ELECTRIC_DIPOLE_DERIVATIVES_MSG << std::endl;
         }
      }
      infile.close();
      return(dmu);
   }
   else
   {
      std::cout << ELECTRIC_DIPOLE_DERIVATIVES_ERR_MSG << std::endl;
      Eigen::MatrixXd dmu;
      return(dmu);
   }
}

Eigen::MatrixXd GaussianFCHK::ReadMagneticDipoleDerivatives() const
{
   if (this->has_magnetic_dipole_derivatives)
   {
      std::ifstream infile;
      infile.open(this->filename);
      const int N = this->Number_of_atoms;
      Eigen::MatrixXd dm = Eigen::MatrixXd::Zero(3,3*N);

      std::string line;
      while (! infile.eof())
      {
         std::getline(infile, line);
         if (line.find(MAGNETIC_DIPOLE_DERIVATIVES_KEY) != std::string::npos) 
         {
            std::vector<double> dm_v = readitems<double>(infile, 3*(3*N));
            for (unsigned int i=0; i < 3*N; ++i) 
            {
	            // the factor 2.0 is required to interpret Gaussian AATs 
               // as magnetic dipole derivatives vs. cartesian nuclear velocities 
               // (see paper by Buckingham, published in Faraday Discussions, 1994)
               dm(0,i) = 2.0 * dm_v[3*i + 0];
               dm(1,i) = 2.0 * dm_v[3*i + 1];
               dm(2,i) = 2.0 * dm_v[3*i + 2];
            }
            std::cout << MAGNETIC_DIPOLE_DERIVATIVES_MSG << std::endl;
         }
      }
      infile.close();
      return(dm);
   }
   else
   {
      std::cout << MAGNETIC_DIPOLE_DERIVATIVES_ERR_MSG << std::endl;
      Eigen::MatrixXd dm;
      return(dm);
   }
}


// diagnostics
bool GaussianFCHK::Has_atomic_numbers() const
{
   return(this->has_atomic_numbers);
}
bool GaussianFCHK::Has_atomic_masses() const
{
   return(this->has_atomic_masses);
}
bool GaussianFCHK::Has_hessian() const
{
   return(this->has_hessian);
}
bool GaussianFCHK::Has_cartesian_coordinates() const
{
   return(this->has_cartesian_coordinates);
}
bool GaussianFCHK::Has_electric_dipole_derivatives() const
{
   return(this->has_electric_dipole_derivatives);
}
bool GaussianFCHK::Has_magnetic_dipole_derivatives() const
{
   return(this->has_magnetic_dipole_derivatives);
}

