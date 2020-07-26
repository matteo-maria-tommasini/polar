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

#ifndef POLARIZABILITY_H
#define POLARIZABILITY_H

#include <complex>
#include <cmath>
#include <Eigen/Core>
#include "states.hpp"

class Polarizability
{

private:
   Eigen::Matrix3cd alpha_ee = Eigen::Matrix3cd::Zero(3,3),
                    alpha_em = Eigen::Matrix3cd::Zero(3,3),
                    alpha_me = Eigen::Matrix3cd::Zero(3,3),
                    alpha_mm = Eigen::Matrix3cd::Zero(3,3);

   double omega = 0.0,
          gamma = 1.0;

   // support function implementing the polarizability formula 
   // (just the contribtion from a given excited state) 
   std::complex<double> alpha_ij_function(const double&, const double&, 
                                          const double&, const double&, 
                                          const double&, const bool&, const bool&) const;

public:
   Polarizability(const States&, const double&, const double&);
   Eigen::Matrix3cd get_ee() const;
   Eigen::Matrix3cd get_mm() const;
   Eigen::Matrix3cd get_me() const;
   Eigen::Matrix3cd get_em() const;
   void Print() const;

};

#endif // POLARIZABILITY_H

