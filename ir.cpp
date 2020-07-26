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

#include "ir.hpp"

IR::IR() : intensities(), dmu_dq()
{
}

IR::IR(const VibraSolver& vibra, const GaussianDataset& gau)
{
   ComputeDipoleDerivativesVsNormalCoordinates(vibra, gau);

   int Nat = gau.GetNumberOfAtoms(); 
   this->intensities = Eigen::MatrixXd::Zero(3*Nat,1);
   for (unsigned int id_mode=0; id_mode < 3 * Nat; ++id_mode)
   {
      // dmu_dq is computed and stored in a.u. (e / A amu^1/2)
      // the factor 2.541765 / A0_RADIUS converts it in debye / (A amu^1/2)
      // IR units: km/mol
      this->intensities(id_mode) = 
                   42.2547 * (
                                pow((2.541765 / A0_RADIUS) * this->dmu_dq[id_mode](0),2.0) + 
                                pow((2.541765 / A0_RADIUS) * this->dmu_dq[id_mode](1),2.0) + 
                                pow((2.541765 / A0_RADIUS) * this->dmu_dq[id_mode](2),2.0) 
                             );
   }
}

void IR::ComputeDipoleDerivativesVsNormalCoordinates(const VibraSolver& vibra, const GaussianDataset& gau)
{
   int Nat = gau.GetNumberOfAtoms(); 

   Eigen::MatrixXd  Lx = vibra.GetCartesianNuclearDisplacements();
   Eigen::MatrixXd dmu =   gau.GetElectricDipoleDerivatives();

   for (unsigned int id_mode = 0; id_mode < 3*Nat; ++id_mode)
   {
      Eigen::Vector3d dmudq = {0.0, 0.0, 0.0}; 
      for (unsigned int id_coord=0; id_coord < 3*Nat; ++id_coord)
      {
          dmudq(0) += dmu(0,id_coord) * Lx(id_coord,id_mode);
          dmudq(1) += dmu(1,id_coord) * Lx(id_coord,id_mode);
          dmudq(2) += dmu(2,id_coord) * Lx(id_coord,id_mode);
      }
      this->dmu_dq.push_back(dmudq);
   }
}

Eigen::MatrixXd IR::GetIntensities() const
{
   return(this->intensities);
}
std::vector<Eigen::Vector3d> IR::GetDipoleDerivativesVsNormalCoordinates() const
{
   return(this->dmu_dq);
}
