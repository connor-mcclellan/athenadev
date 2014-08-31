//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// Primary header
#include "eos.hpp"

// Athena headers
#include "../../athena.hpp"         // enums, macros, Real
#include "../../athena_arrays.hpp"  // AthenaArray
#include "../../mesh.hpp"           // MeshBlock
#include "../fluid.hpp"
#include "../../parameter_input.hpp"

//======================================================================================
/*! \file adiabatic_hydro.cpp
 *  \brief implements functions in class EqnOfState for  adiabatic hydrodynamics`
 *====================================================================================*/

// FluidEqnOfState constructor

FluidEqnOfState::FluidEqnOfState(Fluid *pf, ParameterInput *pin)
{
  pmy_fluid_ = pf;
  gamma_ = pin->GetReal("fluid","gamma");
}

// destructor

FluidEqnOfState::~FluidEqnOfState()
{
}

//--------------------------------------------------------------------------------------
/* \!fn void FluidEqnOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
 *   AthenaArray<Real> &prim_old, AthenaArray<Real> &prim)
 * \brief convert conserved to primitive variables for adiabatic hydro */

void FluidEqnOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
  AthenaArray<Real> &prim_old, AthenaArray<Real> &prim)
{
  MeshBlock *pmb = pmy_fluid_->pmy_block;
  int jl = pmb->js; int ju = pmb->je;
  int kl = pmb->ks; int ku = pmb->ke;
  if (pmb->block_size.nx2 > 1) {
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  if (pmb->block_size.nx3 > 1) {
    kl -= (NGHOST);
    ku += (NGHOST);
  }

  AthenaArray<Real> lcons = cons.ShallowCopy();
  AthenaArray<Real> lprim = prim.ShallowCopy();

//--------------------------------------------------------------------------------------
// Convert to Primitives

  for (int k=kl; k<=ku; ++k){
  for (int j=jl; j<=ju; ++j){
#pragma simd
    for (int i=pmb->is-(NGHOST); i<=pmb->ie+(NGHOST); ++i){
      Real& u_d  = lcons(IDN,k,j,i);
      Real& u_m1 = lcons(IVX,k,j,i);
      Real& u_m2 = lcons(IVY,k,j,i);
      Real& u_m3 = lcons(IVZ,k,j,i);
      Real& u_e  = lcons(IEN,k,j,i);

      Real& w_d  = lprim(IDN,k,j,i);
      Real& w_m1 = lprim(IVX,k,j,i);
      Real& w_m2 = lprim(IVY,k,j,i);
      Real& w_m3 = lprim(IVZ,k,j,i);
      Real& w_p  = lprim(IEN,k,j,i);

      Real di = 1.0/u_d;
      w_d  = u_d;
      w_m1 = u_m1*di;
      w_m2 = u_m2*di;
      w_m3 = u_m3*di;

      w_p = u_e - 0.5*di*(u_m1*u_m1 + u_m2*u_m2 + u_m3*u_m3);
      w_p *= (GetGamma() - 1.0);
    }
  }}

  return;
}