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

#include "states.hpp"
#include "units.hpp"

States::States() : energy(), el_dipole(), mag_dipole()
{
}

unsigned int States::GetNumberOfStates() const
{
   return(this->energy.size());
}

std::vector<Eigen::Vector3d> States::GetElectricDipole() const
{
   return(this->el_dipole);
}

std::vector<Eigen::Vector3d> States::GetMagneticDipole() const
{
   return(this->mag_dipole);
}

void States::DestroyLastNstates(const int N)
{
   int N0 = this->el_dipole.size();
   Eigen::VectorXd resized_energy = this->energy.head(N0-N);
   this->energy = resized_energy;
   for (int i = 0; i<N; ++i)
   {
     this->el_dipole.pop_back();
     this->mag_dipole.pop_back();
   }
}

void States::Print() const
{
   std::cout << "# List of excited states: energies (hartree, eV); transition electric dipole (au); transition magnetic dipole (au):" << std::endl;
   for (int i=0; i<this->energy.size(); ++i) 
   { 
      std::cout << "# " << i+1 << " " 
                << this->energy(i) 
                << " ("
                << this->energy(i) * EV_PER_AU 
                << " eV)"
                << "; " 
                << this->el_dipole[i](0) << " " 
                << this->el_dipole[i](1) << " " 
                << this->el_dipole[i](2) 
                << "; " 
                << this->mag_dipole[i](0) << " " 
                << this->mag_dipole[i](1) << " " 
                << this->mag_dipole[i](2) 

                << std::endl; 
   }

   std::cout << "# List of excited states: energies (eV); oscillatory strength (f); rotatory strength (R)" << std::endl;
   Eigen::VectorXd f = this->OscillatorStrengths(),
                   R = this->RotatoryStrengths();
   for (int i=0; i<this->energy.size(); ++i) 
   { 
      std::cout << "# " << i+1 << " " 
                << this->energy(i) * EV_PER_AU 
                << "; " 
                << f(i) << " " 
                << "; " 
                << R(i) << std::endl; 
   }
}

States::States(const std::string& filename)
{
    std::ifstream infile;
    std::string line;

    int number_of_states = 0;
    int number_of_energy_read = 0;

    infile.open(filename);

    if (infile.fail())
    {
       std::cout << "Unreadable file: " << filename << std::endl;
       std::cout << "Abort." << std::endl;
       exit(0);
    }

    while (! infile.eof())
    {
        std::getline(infile, line);

        // electric transition dipoles
        if (line.find("Ground to excited state transition electric dipole moments (Au):") != std::string::npos) 
        {
           std::getline(infile, line); // this line is skipped
           std::getline(infile, line); 
           while (line.find("Ground to excited state transition velocity dipole moments (Au):") == std::string::npos)
           {
              std::string idx_state; // a string allows dealing properly with Fortran *** overflows
              double mu_x, mu_y, mu_z;
              std::stringstream ss(line);
              ss >> idx_state >> mu_x >> mu_y >> mu_z;
              this->el_dipole.push_back(Eigen::Vector3d(mu_x,mu_y,mu_z));
              // next line
              std::getline(infile, line);
           }
           number_of_states = this->el_dipole.size();
        }
        

        // magnetic transition dipoles
        if (line.find("Ground to excited state transition magnetic dipole moments (Au):") != std::string::npos)
        {
           std::getline(infile, line); // this line is skipped
           for (int i=0; i<number_of_states; ++i)
           {
              std::getline(infile, line);

              std::string idx_state;  // a string allows dealing properly with Fortran *** overflows
              double m_x, m_y, m_z;
              std::stringstream ss(line);
              ss >> idx_state >> m_x >> m_y >> m_z;
              this->mag_dipole.push_back(Eigen::Vector3d(m_x,m_y,m_z));
           }
        }

        // transition energies (converted to hartree)
        if (line.find("Excited State") != std::string::npos)
        {
           std::string         s1, s2, s3, s4;
           double energy; 
           std::stringstream ss(line);
           ss >> s1 >> s2 >> s3 >> s4 >> energy; 
           if (number_of_energy_read == 0) 
           {
              this->energy = Eigen::VectorXd::Zero(number_of_states);  
           }
           this->energy(number_of_energy_read) = energy / EV_PER_AU;
           ++number_of_energy_read;
        }

    }
}


double States::SumRule() const
{
   double sum = 0.0;

   for (int i=0; i<this->energy.size(); ++i)
   {
      // oscillator strength
      double f_i = (2.0/3.0) * this->energy[i] * this->el_dipole[i].squaredNorm(); 
      sum += f_i;
   }
   return(sum);
}

Eigen::VectorXd States::OscillatorStrengths() const
{
   int N = energy.size();
   Eigen::VectorXd f = Eigen::VectorXd::Zero(N);
   for (int i=0; i<this->energy.size(); ++i)
   {
      // oscillator strength
      f(i) = (2.0/3.0) * this->energy[i] * this->el_dipole[i].squaredNorm(); 
   }
   return(f);
}

Eigen::VectorXd States::RotatoryStrengths() const
{
   int N = energy.size();
   Eigen::VectorXd R = Eigen::VectorXd::Zero(N);
   for (int i=0; i<this->energy.size(); ++i)
   {
      // rotatory strength
      R(i) = this->mag_dipole[i].dot(this->el_dipole[i]); 
   }
   return(R);
}

Eigen::VectorXd States::GetElectricTransitionDipoleNorm() const
{
   int N = energy.size();
   Eigen::VectorXd mu = Eigen::VectorXd::Zero(N);
   for (int i=0; i<this->energy.size(); ++i)
   {
      mu(i) = this->el_dipole[i].norm();
   }
   return(mu);
}

Eigen::VectorXd States::GetMagneticTransitionDipoleRobustNorm() const
{
   int N = energy.size();
   Eigen::VectorXd m = Eigen::VectorXd::Zero(N);
   for (int i=0; i<this->energy.size(); ++i)
   {
      // unit vector along the electric transition dipole
      Eigen::Vector3d u = this->el_dipole[i] / this->el_dipole[i].norm();
      // magnetic dipole at the robust point
      m(i) = this->mag_dipole[i].dot(u);
   }
   return(m);
}

Eigen::VectorXd States::GetEnergy() const
{
   return(this->energy);
}
