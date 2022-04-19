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

#include <iostream>
#include <iomanip>
#include <fstream>

#include "icm.hpp"

void ICM::PrintReport() const
{
   if (this->has_ICM_VCD)
   {
      int idx = 1;
      double sum = 0.0;
      for (int i = 0; i<this->ICM_intensity.size(); ++i)
      {
         if (abs(ICM_intensity(i)) > 1.0e-6) 
         {
            std::cout << "VCD lambda " << idx << " = " << ICM_intensity(i) << std::endl;
            ++idx;
            sum += ICM_intensity(i);
         }
      }
      std::cout << "VCD sum of lambdas = " << sum << std::endl; 
      double trace = 0.0;
      for (int i = 0; i<this->ICM_intensity.size(); ++i)
      {
         trace += ICM_intensity(i);
      }
      std::cout << "M-matrix trace = " << trace << std::endl;
   }
   else if (this->has_ICM_IR)
   {
      int idx = 1;
      double sum = 0.0;
      for (unsigned int i = 0; i<this->ICM_intensity.size(); ++i)
      {
         if (abs(ICM_intensity(i)) > 1.0e-6) 
         {
            std::cout << "IR lambda " << idx << " = " << ICM_intensity(i) << std::endl;
            ++idx;
            sum += ICM_intensity(i);
         }
      }
      std::cout << "IR sum of lambdas = " << sum << std::endl; 
      double trace = 0.0;
      for (int i = 0; i<this->ICM_intensity.size(); ++i)
      {
         trace += ICM_intensity(i);
      }
      std::cout << "M-matrix trace = " << trace << std::endl;
   }
   else if (this->has_ICM_MAG)
   {
      int idx = 1;
      double sum = 0.0;
      for (int i = 0; i<this->ICM_intensity.size(); ++i)
      {
         if (abs(ICM_intensity(i)) > 1.0e-6) 
         {
            std::cout << "MAG lambda " << idx << " = " << ICM_intensity(i) << std::endl;
            ++idx;
            sum += ICM_intensity(i);
         }
      }
      std::cout << "MAG sum of lambdas = " << sum << std::endl; 
      double trace = 0.0;
      for (int i = 0; i<this->ICM_intensity.size(); ++i)
      {
         trace += ICM_intensity(i);
      }
      std::cout << "M-matrix trace = " << trace << std::endl;
   }
}

ICM::ICM(const GaussianDataset& gau, const std::string& kind_of_calculation)
{
   // defaults
   this->has_ICM_IR     = false;
   this->has_ICM_VCD    = false;

   // Torii's M-matrix
   Eigen::MatrixXd M;

   /////////////
   // IR ICMs //
   /////////////
   if (kind_of_calculation.compare("IR") == 0)
   {
      if (gau.Has_electric_dipole_derivatives())
      {
         Eigen::MatrixXd dmu = gau.GetElectricDipoleDerivatives();
         int Nat3 = dmu.cols();
         // allocate and form the M-matrix
         M = Eigen::MatrixXd::Zero(Nat3,Nat3);
         for (unsigned int i = 0; i < Nat3; ++i)
         {
            for (unsigned int j = 0; j < Nat3; ++j)
            {
               for (unsigned int k = 0; k < 3; ++k)
               {
                  M(i,j) += dmu(k,i) * dmu(k,j);
               }
            }
         }
         // the end
         this->has_ICM_IR = true;
      }
   }
   //////////////
   // VCD ICMs //
   //////////////
   else if (kind_of_calculation.compare("VCD") == 0)
   {
      if (
            gau.Has_electric_dipole_derivatives() && 
            gau.Has_magnetic_dipole_derivatives()
         )
      {
         Eigen::MatrixXd dmu = gau.GetElectricDipoleDerivatives();
         Eigen::MatrixXd  dm = gau.GetMagneticDipoleDerivatives();
         int Nat3 = dmu.cols();

         // allocate and form the M-matrix
         M = Eigen::MatrixXd::Zero(Nat3,Nat3);
         for (int i = 0; i < Nat3; ++i)
         {
            for (int j = 0; j < Nat3; ++j)
            {
               for (int k = 0; k < 3; ++k)
               {
                  M(i,j) += 0.5 * (dmu(k,i) * dm(k,j) + dmu(k,j) * dm(k,i));
               }
            }
         }
         // the end
         this->has_ICM_VCD = true;
      }
   }
   //////////////
   // MAG ICMs //
   //////////////
   else if (kind_of_calculation.compare("MAG") == 0)
   {
      if (gau.Has_magnetic_dipole_derivatives())
      {
         Eigen::MatrixXd  dm = gau.GetMagneticDipoleDerivatives();
         int Nat3 = dm.cols();

         // allocate and form the M-matrix
         M = Eigen::MatrixXd::Zero(Nat3,Nat3);
         for (int i = 0; i < Nat3; ++i)
         {
            for (int j = 0; j < Nat3; ++j)
            {
               for (int k = 0; k < 3; ++k)
               {
                  M(i,j) += dm(k,i) * dm(k,j);
               }
            }
         }
         // the end
         this->has_ICM_MAG = true;
      }
   }

   // this code is common to any kind of M matrix formed above
   // since ICM are eigenvectors of M
   Eigen::VectorXd eigenvalues;
   Eigen::MatrixXd eigenvectors;
   this->Diagonalize(M, eigenvalues, eigenvectors);
   this->ICM_intensity = eigenvalues;
   this->ICM_displacement = eigenvectors;
}

void ICM::Diagonalize(const Eigen::MatrixXd& W, 
                            Eigen::VectorXd& eigenvalues, 
                            Eigen::MatrixXd& eigenvectors)
{
   // Use SelfAdjointEigenSolver to get eigen values and eigen vectors
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(W);
   Eigen::RowVectorXd unsorted_eigen_values = eigen_solver.eigenvalues();
   Eigen::MatrixXd unsorted_eigen_vectors = eigen_solver.eigenvectors();

   // Stuff below is done to sort eigen values. This can be done in other ways too.
   std::vector<std::pair<int, int>> eigen_value_index_vector;
   for (int i = 0; i < unsorted_eigen_values.size(); ++i)
   {
       eigen_value_index_vector.push_back( std::make_pair(unsorted_eigen_values[i], i));
   }
   std::sort( std::begin(eigen_value_index_vector), 
                std::end(eigen_value_index_vector), 
                std::greater<std::pair<int,int>>() );

   Eigen::RowVectorXd sorted_eigen_values(unsorted_eigen_values.cols());
   Eigen::MatrixXd sorted_eigen_vectors(unsorted_eigen_vectors.rows(), unsorted_eigen_vectors.cols());
   for (int i = 0; i < unsorted_eigen_values.size(); ++i)
   {
       // can also be eigen_value_index_vector[i].first
       sorted_eigen_values[i] = unsorted_eigen_values[eigen_value_index_vector[i].second];
       sorted_eigen_vectors.col(i) = unsorted_eigen_vectors.col(eigen_value_index_vector[i].second);
   }

   // the end
   eigenvalues = sorted_eigen_values;
   eigenvectors = sorted_eigen_vectors;
}

/////////////////
// output ICMs //
/////////////////
void ICM::WriteMolden(const GaussianDataset& gau, const std::string& filename) const 
{
   // 
   //      [Molden Format]
   //  XXX    [SCFCONV]
   //  XXX    [GEOCONV]
   //      [GEOMETRIES] XYZ
   //  XXX    [FORCES]
   //      [FREQ]
   //      [INT]
   //      [FR-COORD]
   //      [FR-NORM-COORD]
   //      

   // Chemical symbols
      char symbtable[][3] = {
      " H", "He", "Li", "Be", " B", " C", " N", " O", " F", "Ne", 
      "Na", "Mg", "Al", "Si", " P", " S", "Cl", "Ar", " K", "Ca", 
      "Sc", "Ti", " V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", 
      "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", " Y", "Zr", 
      "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", 
      "Sb", "Te", " I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", 
      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", 
      "Lu", "Hf", "Ta", " W", "Re", "Os", "Ir", "Pt", "Au", "Hg", 
      "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", 
      "Pa", " U" };

   std::ofstream molden_file;
   molden_file.open(filename);
   if ( (molden_file.rdstate() & std::ifstream::failbit ) != 0 )
   {
      std::cout << "Error opening " << filename << std::endl;
      std::cout << "Abort." << std::endl;
      exit(0);
   }

   int Nat = gau.GetNumberOfAtoms(); 
   int Nmodes = this->ICM_intensity.size();
   if (Nmodes != 3*Nat)
   {
      std::cout << "Error condition in ICM::WriteMolden" << std::endl;
      std::cout << "Abort." << std::endl;
      exit(0);
   }

   std::cout << "Writing Molden-formatted file " << filename << std::endl;

   molden_file << "[Molden Format]"  << std::endl;
   molden_file << "[GEOMETRIES] XYZ" << std::endl;
   molden_file << Nat << std::endl;
   molden_file << "FINAL HEAT OF FORMATION =     0.000000" << std::endl;

   Eigen::MatrixXd coord = gau.GetCartesianCoordinates(); 
   Eigen::MatrixXi atnum = gau.GetAtomicNumbers(); 
   for (int i=0; i<Nat; ++i)
   {
      // coordinates in angstrom
      std::string symbol = symbtable[(int)atnum(i)-1]; 
      molden_file << std::fixed << std::setfill(' ') << std::setw(3) 
                  << symbol << " " 
                  << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(6) 
                  << A0_RADIUS * coord(0,i)     << " " 
                  << A0_RADIUS * coord(1,i)     << " " 
                  << A0_RADIUS * coord(2,i)     << std::endl;
   }

   molden_file << "[FREQ]" << std::endl;
   for (int i=0; i < Nmodes; ++i)
   {
      molden_file << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(2) 
             << this->ICM_intensity(i) << std::endl;
   }

   molden_file << "[INT]" << std::endl;
   for (int i=0; i < Nmodes; ++i)
   {
      molden_file << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(2) 
             << this->ICM_intensity(i) << std::endl;
   }

   molden_file << "[FR-COORD]" << std::endl;
   for (int i=0; i < Nat; ++i)
   {
      std::string symbol = symbtable[(int)atnum(i)-1]; 
      // coordinates in bohr 
      molden_file << std::fixed << std::setfill(' ') << std::setw(3) 
                  << symbol     << " " 
                  << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(6) 
                  << coord(0,i) << " "
                  << coord(1,i) << " " 
                  << coord(2,i) << std::endl;
   }

   molden_file << "[FR-NORM-COORD]" << std::endl;
   for (int idx=0; idx < Nmodes; ++idx)
   {
      molden_file << "vibration" << " " << idx+1 << std::endl;
      for (int i=0; i < Nat; ++i)
      {
         // nuclear displacements 
         molden_file << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(4) 
                     << this->ICM_displacement(3*i + 0, idx) << " " 
                     << this->ICM_displacement(3*i + 1, idx) << " "
                     << this->ICM_displacement(3*i + 2, idx) << std::endl;
      }
   }

   // close file
   molden_file.close();
}


