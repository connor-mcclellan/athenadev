// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code
// contributors Licensed under the 3-clause BSD License, see LICENSE file for
// details
//========================================================================================
//! \file eddington_wind.cpp
//  \brief 1D wind driven by a constant flux inner boundary condition
//========================================================================================

// C++ headers
#include <algorithm> // min
#include <cfloat>    // FLT_MAX
#include <cmath>     // sqrt
#include <cstdlib>   // srand
#include <fstream>
#include <iomanip>
#include <iostream>  // endl
#include <sstream>   // stringstream
#include <stdexcept> // runtime_error
#include <string>    // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../inputs/mesa_reader.hpp"
#include "../mesh/mesh.hpp"
#include "../nr_radiation/integrators/rad_integrators.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/utils.hpp"

namespace {

// Dimensional quantities
Real r0;        // Length
Real v0;        // Velocity
Real rho0;      // Density
Real temp0;     // Temperature

// Boundary conditions
Real vb;        // Inflow velocity
Real rhob;      // Base density
Real tempb;     // Base temperature

// Physical parameters
Real mu;        // mean molecular weight
Real gm;        // gravitational parameter

Real mdotb, egamb;

} // namespace

void UserPointMass(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
    AthenaArray<Real> &cons_scalar);

void HydroInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt, int il, int iu, int jl,
                  int ju, int kl, int ku, int ngh);


void Mesh::InitUserMeshData(ParameterInput *pin) {

  v0 = pin->GetReal("problem", "v0");
  temp0 = pin->GetReal("problem", "temp0");
  r0 = pin->GetReal("problem", "r0");
  rho0 = pin->GetReal("problem", "rho0");

  gm = pin->GetReal("problem", "gm");
  mu = pin->GetReal("problem", "mu");
  rhob = pin->GetReal("problem", "rhob");
  vb = pin->GetReal("problem", "vb");
  tempb = pin->GetReal("problem", "tempb");

  // Mass loss rate at inner boundary
  Real rb = pin->GetReal("mesh", "x1min");
  Real gamma = pin->GetReal("hydro", "gamma");
  mdotb = vb * 4.0 * PI * SQR(rb) * rhob;
  egamb = 0.5*SQR(vb) + (gamma / (gamma-1.0)) * tempb - gm / rb;

  EnrollUserExplicitSourceFunction(UserPointMass);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, HydroInnerX1);
  return;
} // InitUserMeshData


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {

  // Get the size of the mesh, including ghost zones
  int nx1 = pmy_mesh->mesh_size.nx1 + 2 * NGHOST;
  int nx2 = pmy_mesh->mesh_size.nx2 + 2 * NGHOST;
  int nx3 = pmy_mesh->mesh_size.nx3 + 2 * NGHOST;

  AllocateRealUserMeshBlockDataField(Uov::NUM_UOV);
  for (int l=0; l<Uov::NUM_UOV; l++) {
    ruser_meshblock_data[l].NewAthenaArray(nx3, nx2, nx1);
  }
  AllocateUserOutputVariables(Uov::NUM_UOV);

  SetUserOutputVariableName(Uov::CSOUND, "csound");
  SetUserOutputVariableName(Uov::GRAD_P, "grad_P");
  SetUserOutputVariableName(Uov::GRAVSRC, "gravsrc");
  SetUserOutputVariableName(Uov::MDOT, "mdot");
  SetUserOutputVariableName(Uov::EGAM, "egam");

} // InitUserMeshblockData


void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {

  AthenaArray<Real> &x1flux = phydro->flux[X1DIR];
  Real gamma = peos->GetGamma();

  // Copy over meshblock data for user output variables
  for (int k = ks; k <= ke; k++) {
    for (int j = js; j <= je; j++) {
      for (int i = is - NGHOST; i <= ie + NGHOST; i++) {
        for (int l = 0; l < Uov::NUM_UOV; l++) {

          if (l==Uov::GRAD_P) {

            if ((i!=0) && (i!=ie+NGHOST)) {
              // Numerical estimate of the change in pressure across this cell
              Real press_l = phydro->w(IPR,k,j,i-1);
              Real press_c = phydro->w(IPR,k,j,i);
              Real press_r = phydro->w(IPR,k,j,i+1);
              Real diff_P = 0.5*(press_r+press_c) - 0.5*(press_c+press_l);
              user_out_var(Uov::GRAD_P,k,j,i) = diff_P / pcoord->dx1f(i);
            } else {
              user_out_var(Uov::GRAD_P,k,j,i) = 0.0;
            }

          } else if (l==Uov::CSOUND) {
            Real rho = phydro->u(IDN, k, j, i);
            Real press = phydro->w(IPR, k, j, i);
            user_out_var(Uov::CSOUND, k, j, i) = std::sqrt(gamma * press / rho);

          } else if (l==Uov::MDOT) {
            Real x1v = pcoord->x1v(i);
            Real rho = phydro->u(IDN, k, j, i);
            Real vel = phydro->w(IVX, k, j, i);
            user_out_var(Uov::MDOT, k, j, i) = (4.0 * PI * SQR(x1v) * rho * vel);

          } else if (l==Uov::EGAM) {
            Real x1v = pcoord->x1v(i);
            Real rho = phydro->u(IDN, k, j, i);
            Real vel = phydro->w(IVX, k, j, i);
            Real press = phydro->w(IPR, k, j, i);
            user_out_var(Uov::EGAM, k, j, i) = (0.5*SQR(vel) + (gamma / (gamma-1.0))
                                                * press / rho - gm / x1v);

          } else {
            user_out_var(l, k, j, i) = ruser_meshblock_data[l](k, j, i);
          }
        }
      }
    }
  }
} // UserWorkBeforeOutput


void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real gamma = peos->GetGamma();
  Real x1min = pmy_mesh->mesh_size.x1min;

  // Initialize the gas
  for (int k = ks; k <= ke; k++) {
    for (int j = js; j <= je; j++) {
      for (int i = is - NGHOST; i <= ie + NGHOST; i++) {

        // Initialize velocity consistent with a constant mass loss rate
        Real x1v = pcoord->x1v(i);
        Real rho = rhob * std::pow(x1v / x1min, -2.0);
        Real vel = mdotb / (4.0 * PI * SQR(x1v) * rho);

        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = vel * rho;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        // Lamers & Cassinelli (1999), Eq. 4.32
        Real temp = tempb * std::pow(rho / rhob, gamma - 1.0);

        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN, k, j, i) = rho * temp / (gamma - 1.0);
          phydro->u(IEN, k, j, i) += 0.5 * SQR(phydro->u(IM1, k, j, i)) / rho;
          phydro->u(IEN, k, j, i) += 0.5 * SQR(phydro->u(IM2, k, j, i)) / rho;
          phydro->u(IEN, k, j, i) += 0.5 * SQR(phydro->u(IM3, k, j, i)) / rho;
        }
      }
    }
  }
}

void UserPointMass(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
    AthenaArray<Real> &cons_scalar) {
  // Copies the built-in PointMass hydro source term exactly for spherical
  // polar coordinates. Used to update the user meshblock data even when
  // hydro evolution is turned off.

  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  Coordinates *pcoord = pmb->pcoord;
  AthenaArray<Real> &x1flux = pmb->phydro->flux[X1DIR];

  for (int i = is; i <= ie; ++i) {
    for (int j = js; j <= je; ++j) {
      for (int k = ks; k <= ke; ++k) {

        // Mass density
        Real rho = prim(IDN, k, j, i);

        // Left and right cell faces
        Real rl = pcoord->x1f(i);
        Real rr = pcoord->x1f(i + 1);
        Real rc = pcoord->x1v(i);

        // Gravitational potential at cell center, left face, and right face
        Real phic = -gm/pcoord->x1v(i);
        Real phil = -gm/rl;
        Real phir = -gm/rr;

        Real vol = (rr*rr*rr - rl*rl*rl) / 3.0;
        Real area = 0.5*(rr*rr - rl*rl);
        Real coord_src = area / vol;

        Real src = dt * rho * coord_src * gm / rc;
        pmb->ruser_meshblock_data[Uov::GRAVSRC](k, j, i) = -src / dt;

        // only update the conserved variables if fluid evolution is turned on
        if (pmb->pmy_mesh->fluid_setup == FluidFormulation::evolve) {
          cons(IM1, k, j, i) -= src;
          if (NON_BAROTROPIC_EOS) {
            Real phy_src = 1.0 / SQR(rc);
            cons(IEN, k, j, i) -= dt * 0.5 * gm * phy_src
                                  * (x1flux(IDN,k,j,i) + x1flux(IDN,k,j,i+1));
          }
        }
      }
    }
  }
}


void HydroInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt, int il, int iu, int jl,
                  int ju, int kl, int ku, int ngh) {

  Real gamma = pmb->peos->GetGamma();
  Real x1min = pmb->pmy_mesh->mesh_size.x1min;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il-1; i>=il-ngh; --i) {
        for (int n=0; n<=NHYDRO; n++) {
          Real x1v = pmb->pcoord->x1v(i);
          if (n == IDN) {
            Real rho = rhob * std::pow(x1v / x1min, -2.0);
            prim(n, k, j, i) = rho;
          } else if (n == IPR) {
            Real rho = prim(IDN,k,j,i);
            prim(n, k, j, i) = rho * tempb * std::pow(rho / rhob, gamma - 1.0);
          } else if (n == IVX) {
            prim(n, k, j, i) = vb;
          } else {
            prim(n, k, j, i) = TINY_NUMBER;
          }
        }
      }
    }
  }
  return;
}


namespace {
} // namespace
