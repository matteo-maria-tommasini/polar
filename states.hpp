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

#ifndef STATES_H
#define STATES_H

#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <vector>

class Polarizability;

class States
{

friend class Polarizability;

private:
   // attributes
   Eigen::VectorXd energy;
   std::vector<Eigen::Vector3d>  el_dipole,
                                mag_dipole;

public:
   States();
   explicit States(const std::string&);
   void Print() const;
   Eigen::VectorXd OscillatorStrengths() const;
   Eigen::VectorXd RotatoryStrengths() const;
   double SumRule() const;
   void DestroyLastNstates(const int);

   // getters
   unsigned int GetNumberOfStates() const;
   std::vector<Eigen::Vector3d> GetElectricDipole() const;
   std::vector<Eigen::Vector3d> GetMagneticDipole() const;
   Eigen::VectorXd GetElectricTransitionDipoleNorm() const;
   Eigen::VectorXd GetMagneticTransitionDipoleRobustNorm() const;
   Eigen::VectorXd GetEnergy() const;

};

#endif // STATES_H

