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

#ifndef IR_H
#define IR_H

#include "units.hpp"
#include "gaussiandataset.hpp"
#include "vibrasolver.hpp"

/////////////////////////////////////////////////////////////////////////////////
// IR is a legacy class meant to reproduce data provided in Gaussian09 output, //
// namely IR absorptions in km/mol. The simulation of IR spectra is taken care //
// of by the vibpolarizability class.                                          //
/////////////////////////////////////////////////////////////////////////////////
class IR 
{
   private:

   Eigen::MatrixXd intensities;
   std::vector<Eigen::Vector3d> dmu_dq;

   public:

   IR();
   IR(const VibraSolver&, const GaussianDataset&);

   // getters
   Eigen::MatrixXd GetIntensities() const;
   std::vector<Eigen::Vector3d> GetDipoleDerivativesVsNormalCoordinates() const;

   // utility
   void ComputeDipoleDerivativesVsNormalCoordinates(const VibraSolver&, const GaussianDataset&);

};

#endif // IR_H

