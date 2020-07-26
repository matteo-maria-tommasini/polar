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

#ifndef ICM_H
#define ICM_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues> 

#include "units.hpp"
#include "options.hpp"
#include "gaussiandataset.hpp"

///////////////////////////////////
// ICM: Intensity Carrying Modes //
///////////////////////////////////
class ICM 
{
   private:

   bool has_ICM_IR  = false,
        has_ICM_VCD = false,
        has_ICM_MAG = false;

   Eigen::VectorXd ICM_intensity;
   Eigen::MatrixXd ICM_displacement;

   public:
   ICM() = default;
   explicit ICM(const GaussianDataset&, const std::string&);
   void Diagonalize(const Eigen::MatrixXd&, Eigen::VectorXd&, Eigen::MatrixXd&);
   void WriteMolden(const GaussianDataset&, const std::string&) const; 
   void PrintReport() const;

};

#endif // ICM_H

