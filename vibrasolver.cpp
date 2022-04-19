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

#include "vibrasolver.hpp"
#include "ir.hpp"
#include "vcd.hpp"
#include "macros.hpp"
#include "units.hpp"
#include "options.hpp"

std::vector<Eigen::Vector3d> VibraSolver::GetElectricDipoleDerivativesVsNormalCoordinates() const
{
   return(this->dmu_dq);
}

std::vector<Eigen::Vector3d> VibraSolver::GetMagneticDipoleDerivativesVsNormalCoordinates() const
{
   return(this->dm_dqdot);
}

void VibraSolver::Diagonalize(const Eigen::MatrixXd& W, 
                              Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors)
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

   eigenvalues = sorted_eigen_values;
   eigenvectors = sorted_eigen_vectors;
}
 

void VibraSolver::SolveSecularEquation(const GaussianDataset& gau, 
                                       Eigen::MatrixXd& Lx, Eigen::VectorXd& wavenumbers)
{
   int Nat = gau.GetNumberOfAtoms();
   Eigen::MatrixXd mass = gau.GetAtomicMasses(); 
   Eigen::MatrixXd invM05 = Eigen::MatrixXd::Zero(3*Nat,3*Nat); 
   Eigen::MatrixXd M = Eigen::MatrixXd::Zero(3*Nat,3*Nat); 

   for (int i = 0; i < Nat; ++i)
   {
      double m05 = sqrt(mass(i));
      invM05(3*i+0,3*i+0) = 1.0 / m05;
      invM05(3*i+1,3*i+1) = 1.0 / m05;
      invM05(3*i+2,3*i+2) = 1.0 / m05;

      M(3*i+0,3*i+0) = mass(i);
      M(3*i+1,3*i+1) = mass(i);
      M(3*i+2,3*i+2) = mass(i);
   }

   // numerical value: const1 = 1302.7909462365
   double const1 = std::sqrt(NAV * 1.0e5) / (2.0 * PI * SPEEDC_CGS); 
   double const2 = 4.3597482 / (std::pow(A0_RADIUS,2)); 

   Eigen::MatrixXd W = const2 * invM05 * gau.GetHessian() * invM05;
   Eigen::VectorXd eigenvalues;
   Eigen::MatrixXd eigenvectors;
   Diagonalize(W, eigenvalues, eigenvectors);

   // wavenumbers: initialization
   wavenumbers = Eigen::VectorXd::Zero(3*Nat); 
   for (int i = 0; i < 3*Nat; ++i)
   {
      wavenumbers(i) = eigenvalues(i);
   }

   // wavenumbers: fix units and negative eigenvalues
   for (int i = 0; i < 3*Nat; ++i)
   {
      if (wavenumbers(i) > 0.0) 
          wavenumbers(i) =  const1 * std::sqrt( wavenumbers(i));
      else
          wavenumbers(i) = -const1 * std::sqrt(-wavenumbers(i));
   }

   // nuclear displacements Lx
   Lx = invM05 * eigenvectors;

}

VibraSolver::VibraSolver(const GaussianDataset& gau)
{
   // defaults
   this->has_IR     = false;
   this->has_VCD    = false;

   SolveSecularEquation(gau, this->Lx, this->wavenumbers); 

   ////////////////////
   // IR intensities //
   ////////////////////
   if (gau.Has_electric_dipole_derivatives())
   {
      IR ir = IR((*this), gau);
      this->IR_intensities = ir.GetIntensities();
      this->dmu_dq = ir.GetDipoleDerivativesVsNormalCoordinates();
      this->has_IR = true;
   }

   /////////////////////
   // VCD intensities //
   /////////////////////
   if (
         gau.Has_electric_dipole_derivatives() && 
         gau.Has_magnetic_dipole_derivatives()
      )
   {
      VCD vcd = VCD((*this), gau);
      this->VCD_intensities = vcd.GetIntensities();
      this->dm_dqdot = vcd.GetMagneticDipoleDerivativesVsNormalCoordinates();
      this->m10_robust = vcd.Get_m10_Robust();
      this->mu01_robust = vcd.Get_mu01_Robust();
      this->has_VCD = true;
   }

}

// getters
int VibraSolver::GetNumberOfAtoms() const
{
   return(this->Lx.rows() / 3);
}
Eigen::MatrixXd VibraSolver::GetCartesianNuclearDisplacements() const
{
   return(this->Lx);
}
Eigen::VectorXd VibraSolver::GetWavenumbers() const
{
   return(this->wavenumbers);
}
Eigen::VectorXd VibraSolver::GetIRIntensities() const
{
   return(this->IR_intensities);
}
Eigen::VectorXd VibraSolver::GetVCDIntensities() const
{
   return(this->VCD_intensities);
}

Eigen::VectorXd VibraSolver::Get_m10_Robust() const
{
   return(this->m10_robust);
}

Eigen::VectorXd VibraSolver::Get_mu01_Robust() const
{
   return(this->mu01_robust);
}

////////////////////////
// main output method //
////////////////////////
void VibraSolver::PrintReport() const
{
   ///////////////////
   // HEADER begins //
   ///////////////////
   std::cout << BAR << std::endl;
   std::cout << "Number of Cartesian degrees of freedom (3N): " 
             << this->wavenumbers.size() << std::endl;
   std::cout << "3N-6 = " 
             << this->wavenumbers.size() - 6 << std::endl;
   std::cout << BAR << std::endl;

   std::vector<std::string> labels;
                         labels.push_back("wavenumber, 1/cm");
   if (this->has_IR)     labels.push_back("IR intensity, km/mol");
   if (this->has_VCD)    labels.push_back("VCD intensity (rotatory strength), 1.0e-44 esu**2 cm**2");
   
   {
      auto N = labels.size(), i = N;
      for (i=0; i<labels.size(); ++i)
      {
         std::cout << "(" << i+1 << ") " << labels[i] << std::endl; 
      }
      std::cout << BAR << std::endl;

      for (i=0; i<labels.size(); ++i)
      {
         std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(2)
                   << i+1;
      }
   }
   std::cout << std::endl;
   std::cout << BAR << std::endl;
   /////////////////
   // HEADER ends //
   /////////////////
  
   ///////////////////////////////////
   // print data rows, mode by mode //
   ///////////////////////////////////
   {
      std::vector<double> values;
      auto N3 = this->wavenumbers.size(), i = N3;
      auto Nval = values.size(), k = Nval;
      for (i=0; i<N3; ++i)
      {
                               values.push_back(this->wavenumbers(i));
         if (this->has_IR)     values.push_back(this->IR_intensities(i));
         if (this->has_VCD)    values.push_back(this->VCD_intensities(i) * 1.0e44);
         
         // final printout
         for (k=0; k<Nval; ++k)
         {
            std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(2)
                      << values[k];
         }
         std::cout 
                   << std::fixed << std::setfill(' ') << std::setw(6) 
                   << "MODE" 
                   << std::fixed << std::setfill(' ') << std::setw(5) 
                   << (i+1)                           
                   << std::endl;

         values.clear();
      }
   }

}

