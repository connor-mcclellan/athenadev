// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file eddington_wind.cpp
//  \brief 1D wind driven by a constant flux inner boundary condition
//========================================================================================

// C++ headers
#include <algorithm>  // min
#include <cfloat>     // FLT_MAX
#include <cmath>      // sqrt
#include <cstdlib>    // srand
#include <fstream>
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../nr_radiation/integrators/rad_integrators.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../utils/utils.hpp"

namespace {

  // Global variables (scope includes this file only)
  // Note that the BC functions **do not** have access to these

  Real dens_base; // Density at x1min (will be higher in ghost zones)
  Real vinflow;   // Velocity initialized everywhere and continually enforced in inner
                  // ghost zones

  // CHARACTERISTIC DIMENSIONAL QUANTITIES

  Real mdot;      // Mass loss rate
  Real edot;      // Energy loss rate
  Real ledd;      // Eddington luminosity
  Real kappa0;    // Opacity
  Real temp0;     // Temperature
  Real rho0;      // Density
  Real r0;        // Length

  // Coordinate information for user output variables
  AthenaArray<Real> x1area;
  AthenaArray<Real> vol;

  // Utility function for getting lab frame intensity from comoving energy and flux
  void SetLabIr(Real *lab_ir, Real *mu, Real *wmu, int nang, Real Er_com, Real Fr_com,
      Real beta, int k, int j, int i);
}

// Function Declarations
void UserOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);

void Newtonian(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
    AthenaArray<Real> &cons_scalar);

void RadFixedInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
    const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir, Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void RadVacuumOuterX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
    const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir, Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void FixedInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void VacuumOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh);


//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Physical constants
  Real c_cgs = 2.997924589e10;
  Real kb    = 1.380649e-16;
  Real mp    = 1.6726e-24;
  Real GMsun = 1.32712440018e26;

  // Simulation input parameters
  Real crat = pin->GetReal("radiation", "crat");
  Real prat = pin->GetReal("radiation", "prat");

  // Problem-specific constants
  Real mmw = pin->GetReal("problem", "mean_mol_weight");
  Real v0 = c_cgs / crat;
  Real kappa_es = 0.2;

  // File-scope problem input parameters
  dens_base = pin->GetReal("problem", "base_density");

  // Characteristic dimensional quantities
  Real mass = pin->GetReal("problem", "mass");
  r0 = mass * GMsun / SQR(c_cgs);
  kappa0 = kappa_es;
  rho0 = 1. / kappa0 / r0;
  temp0 = SQR(v0) * mmw * mp / kb;

  mdot = pin->GetReal("problem", "mdot") / (rho0*v0*SQR(r0));
  edot = pin->GetReal("problem", "edot");

  ledd = 2*PI*pow(crat,3);
  vinflow = mdot/(4.*PI*dens_base*SQR(mesh_size.x1min));

  // Enroll boundary functions and gravitational source term
  EnrollUserExplicitSourceFunction(Newtonian);
  EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, RadFixedInnerX1);
  EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, RadVacuumOuterX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, FixedInnerX1);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, VacuumOuterX1);

  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {

  pnrrad->EnrollOpacityFunction(UserOpacity);

  int nx1 = pmy_mesh->mesh_size.nx1+2*NGHOST;
  int nx2 = pmy_mesh->mesh_size.nx2+2*NGHOST;
  int nx3 = pmy_mesh->mesh_size.nx3+2*NGHOST;

  x1area.NewAthenaArray(nx1+1);
  vol.NewAthenaArray(nx1);
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pcoord->Face1Area(k, j, is-NGHOST, ie+NGHOST+1, x1area);
      pcoord->CellVolume(k, j, is-NGHOST, ie+NGHOST, vol);
    }
  }

  // Allocate user meshblock data arrays
  AllocateRealUserMeshBlockDataField(2);
  ruser_meshblock_data[0].NewAthenaArray(nx3, nx2, nx1);
  ruser_meshblock_data[1].NewAthenaArray(nx3, nx2, nx1);

  // Allocate user output variables and set their names
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "GravSrc_IM1");
  SetUserOutputVariableName(1, "RadSource1");
  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {

  Real gamma = peos->GetGamma();
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pcoord->x1f(ie+NGHOST+1);

  for (int k = ks; k <= ke; k++) {
    for (int j = js; j <= je; j++) {
      for (int i = is-NGHOST; i <= ie+NGHOST; i++) {

        Real rho = phydro->u(IDN, k, j, i);

        // Calculate gravity source term

        // Left and right cell faces
        Real rl = pcoord->x1f(i);
        Real rr = pcoord->x1f(i + 1);
        Real rc = pcoord->x1v(i);

        // Gravitational parameter
        Real gm = 0.5 * SQR(pnrrad->crat);

        // Gravitational potential
        Real phil = -gm/rl;
        Real phir = -gm/rr;

        // Source term - force per unit volume
        Real dt = pmy_mesh->dt;
        Real src = - dt * rho * (phir - phil) / (rr - rl);

        // Gas acceleration per timestep due to gravity
        user_out_var(0, k, j, i) = src / dt / rho;

        // Copy over meshblock data for RadSource1
        user_out_var(1, k, j, i) = ruser_meshblock_data[0](k, j, i) / dt / rho;
      }
    }
  }
}
//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the X-ray bursts problem
//========================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pcoord->x1f(ie+NGHOST+1); // TODO: This could fail if nblocks > 1?

  Real Fr_com_target = 0.5 * SQR(pnrrad->crat) / pnrrad->prat / (SQR(x1min));

  for (int k = ks; k <= ke; k++) {
    for (int j = js; j <= je; j++) {
      for (int i = is-NGHOST; i <= ie+NGHOST; i++) {

        // Initialize powerlaw density profile
        Real x = pcoord->x1v(i);
        phydro->u(IDN, k, j, i) = dens_base * pow(x/x1min, -2);

        // Initialize momentum
        phydro->u(IM1, k, j, i) = phydro->u(IDN,k,j,i) * vinflow;
        phydro->u(IM2, k, j, i) = 0.0;
        phydro->u(IM3, k, j, i) = 0.0;
      }
    }
  }

  if (NON_BAROTROPIC_EOS) {
    Real gamma = peos->GetGamma();
    for (int k = ks; k <= ke; k++) {
      for (int j = js; j <= je; j++) {
        for (int i = is-NGHOST; i <= ie+NGHOST; i++) {

          Real rho = phydro->u(IDN, k, j, i);

          // Get gas temperature that will put the gas and radiation
          // roughly in equilibrium
          Real x = pcoord->x1v(i);
          Real tau = SQR(x1min) * dens_base * (1./x - 1./x1max);
          Real Er_com_target = 3. * Fr_com_target * tau;
          Real temp = std::pow(Er_com_target, 0.25);

          // Initialize gas internal energy
          phydro->u(IEN, k, j, i) = rho * temp / (gamma - 1.0);

          // Add on gas kinetic energy
          phydro->u(IEN, k, j, i) += 0.5 * SQR(phydro->u(IM1, k, j, i)) / rho;
          phydro->u(IEN, k, j, i) += 0.5 * SQR(phydro->u(IM2, k, j, i)) / rho;
          phydro->u(IEN, k, j, i) += 0.5 * SQR(phydro->u(IM3, k, j, i)) / rho;
        }
      }
    }
  }

  if (IM_RADIATION_ENABLED) {
    for (int k = ks; k <= ke; k++) {
      for (int j = js; j <= je; j++) {
        for (int i = is-NGHOST; i <= ie+NGHOST; i++) {

          // Set intensities everywhere to yield a comoving radial flux equal
          // to Fr_com_target, with a consistent energy density under the Eddington
          // approximation (Pr = 1/3 Er) and at large optical depth (c Er = 3 Fr \tau)

          Real x = pcoord->x1v(i);
          Real tau = SQR(x1min) * dens_base * (1./x - 1./x1max);
          Real Er_com_target = 3. * Fr_com_target * tau;

          for (int ifr = 0; ifr < pnrrad->nfreq; ++ifr) {
            int n0 = ifr*pnrrad->nang;
            Real *lab_ir = &(pnrrad->ir(k, j, i, n0));
            Real *mu     = &(pnrrad->mu(0, k, j, i, n0));
            Real *wmu    = &(pnrrad->wmu(n0));
            SetLabIr(lab_ir, mu, wmu, pnrrad->nang, Er_com_target, Fr_com_target,
                vinflow/pnrrad->crat, k, j, i);
          }
        }
      }
    }
  }
  return;
}


void Newtonian(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
    AthenaArray<Real> &cons_scalar) {

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

        // Gravitational parameter
        Real gm = 0.5 * SQR(pmb->pnrrad->crat);

        // Gravitational potential at cell center, left face, and right face
        Real phic = -gm/pcoord->x1v(i);
        Real phil = -gm/rl;
        Real phir = -gm/rr;

        // Cell-average position cubed
        Real vol = (rr*rr*rr - rl*rl*rl)/3;

        // Source term - force per unit volume
        Real src = - dt * rho * (phir - phil) / (rr - rl);

        cons(IM1, k, j, i) += src;
        pmb->ruser_meshblock_data[0](k, j, i) = src / dt / rho;

        // 1D cell face areas
        Real arear = rr * rr;
        Real areal = rl * rl;

        // x1flux is mass flux - g / (cm^2 s)
        // area * x1flux is mass loss rate - g / s
        Real phidivrhov = (arear * x1flux(IDN, k, j, i+1) -
                           areal * x1flux(IDN, k, j, i)) * phic/vol;

        Real divrhovphi = (arear * x1flux(IDN, k, j, i+1) * phir -
                           areal * x1flux(IDN, k, j, i) * phil) / vol;

        cons(IEN, k, j, i) += dt * (phidivrhov - divrhovphi);
      }
    }
  }
}

void UserOpacity(MeshBlock *pmb, AthenaArray<Real> &prim) {

  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  NRRadiation *pnrrad = pmb->pnrrad;

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is-NGHOST; i <= ie+NGHOST; ++i) {
        for (int ifr = 0; ifr < pnrrad->nfreq; ++ifr) {

          Real rho = prim(IDN,k,j,i);
          Real temp = prim(IPR,k,j,i) / rho;

          // Free-free opacity from (???) (shane calculated this)
          Real kappa_ff = 2.87e24 * (rho * rho0) / pow(temp * temp0, 3.5); // Check R&L
          kappa_ff /= kappa0; // dimensionless

          // Scattering opacity
          Real kappa_s =  1.; // dimensionless

          // Set all the matter-radiation coupling coefficients
          pnrrad->sigma_s(k, j, i, ifr) = rho * kappa_s;
          pnrrad->sigma_a(k, j, i, ifr) = 0.0;
          pnrrad->sigma_p(k, j, i, ifr) = rho * kappa_ff;
          pnrrad->sigma_pe(k, j, i, ifr) = rho * kappa_ff;
        }
      }
    }
  }
}

// Constant flux through inner x1 boundary
void RadFixedInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
    const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir, Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  Real gamma = pmb->peos->GetGamma();
  Real gfac = gamma/(gamma-1);

  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {

      // Get the lab frame energy density in first active zone
      Real Er = 0.0;
      for (int ifr=0; ifr<pnrrad->nfreq; ++ifr) {
        Real *lab_ir = &(ir(k, j, il, ifr*pnrrad->nang));
        Real *wmu    = &(pnrrad->wmu(ifr*pnrrad->nang));

        for (int n=0; n<pnrrad->nang; n++) {
          Er += lab_ir[n] * wmu[n];
        }
      }

      // Set intensity in ghost zones
      for (int i=0; i<NGHOST; ++i) {
	Real x1v = pmb->pcoord->x1v(il-1-i);
	Real Fr_edd = 0.5 * SQR(pnrrad->crat) / pnrrad->prat / (SQR(x1v));
	Real cs2 = pow(Er,0.25)/w(IDN,k,j,i);
	Real bern = 0.5*(SQR(vinflow)-SQR(pnrrad->crat)/x1v)+gfac*cs2;
	Real Fr_target = Fr_edd*(edot-mdot/(2.*PI*pow(pnrrad->crat,3))*bern);
        for (int ifr=0; ifr<pnrrad->nfreq; ++ifr) {
          Real *lab_ir = &(ir(k, j, i, ifr*pnrrad->nang));
          Real *mu     = &(pnrrad->mu(0, k, j, i, ifr*pnrrad->nang));
          Real *wmu    = &(pnrrad->wmu(ifr*pnrrad->nang));
          SetLabIr(lab_ir, mu, wmu, pnrrad->nang, Er, Fr_target, 0., k, j, i);
        }
      }
    }
  }
}

void RadVacuumOuterX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
    const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir, Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  // copy radiation variables into ghost zones
  int nang = pnrrad->nang; // number angles
  int nfreq = pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for (int n=0; n<nang; ++n) {
            int ang=ifr*nang+n;
	    Real mux= pnrrad->mu(0,k,j,iu,n);
	    if(mux > 0.0){
	      ir(k,j,iu+i,ang) = ir(k,j,iu,ang);
	    }else{
	      ir(k,j,iu+i,ang) = 0.0;
	    }
          }
        }
      }
    }
  }
  return;
}

void FixedInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
              Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku,
              int ngh) {

  NRRadiation *pnrrad = pmb->pnrrad;
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      Real Er = 0.0;
      for (int ifr=0; ifr<pnrrad->nfreq; ++ifr) {
        Real *lab_ir = &(pnrrad->ir(k, j, il, ifr*pnrrad->nang));
        Real *wmu    = &(pnrrad->wmu(ifr*pnrrad->nang));

        for (int n=0; n<pnrrad->nang; n++) {
          Er += lab_ir[n] * wmu[n];
        }
      }
      for (int i=1; i<=ngh; ++i) {
        for (int n=0; n<=NHYDRO; n++) {
          if (n == IDN) {

              prim(n, k, j, il-i) = dens_base;

            } else if (n == IPR) {

              Real x = pmb->pcoord->x1v(i);
	      Real tgas = pow(Er,0.25);
              prim(n, k, j, il-i) = prim(IDN, k, j, il) * tgas;

            } else if (n == IVX) {
              prim(IVX, k, j, il-i) = vinflow;
            } else {
              prim(n, k, j, il-i) = prim(n, k, j, il);
            }
          }
        }
      }
  }
  return;
}

void VacuumOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
  FaceField &b, Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku,
  int ngh) {

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for (int n=0; n<=NHYDRO; n++) {
          if (prim(IVX,k,j,iu) >= 0.) {
            prim(n,k,j,iu+i) = prim(n,k,j,iu);
          } else {
            prim(n,k,j,iu+i) = 0.0;
          }
        }
      }
    }
  }
  return;
}





namespace {

  Real LabOverCom(Real mu, Real beta) {
    Real lorz = sqrt(1.0 / (1.0 - SQR(beta)));
    Real bdotn = beta * mu;
    return SQR(lorz * (1.0 - bdotn));
  }

  // Set lab frame intensity using a target energy density and flux in the comoving frame
  void SetLabIr(Real *lab_ir, Real *mu, Real *wmu, int nang, Real Er_com, Real Fr_com,
      Real beta, int k, int j, int i) {

    Real er_pos, fr_pos;
    Real er_neg, fr_neg;
    Real wnorm_com;

    // Ensure that weights sum to 1 in the comoving frame
    wnorm_com = 0.;
    for (int n=0; n<nang; ++n) {
      wnorm_com += wmu[n] / LabOverCom(mu[n], beta);
    }

    // Sum upward and downward intensity coefficients
    er_pos = 0., fr_pos = 0.;
    er_neg = 0., fr_neg = 0.;
    for (int n=0; n<nang; ++n) {
      Real lab_over_com = LabOverCom(mu[n], beta);

      // Weights in comoving frame
      Real wmu_com = wmu[n] / lab_over_com / wnorm_com;

      // Transform angle to comoving frame
      // Mihalas & Mihalas, Foundations of Radiation Hydrodynamics, Eq 89.6
      Real lorz = sqrt(1.0 / (1.0 - SQR(beta)));
      Real mu_com = (mu[n] - (lorz * beta) * (1.0 - lorz / (lorz + 1.0)
                    * (beta * mu[n]))) / sqrt(lab_over_com);

      if (mu_com > 0.0) {
        er_pos += wmu_com;
        fr_pos += (wmu_com * mu_com);
      } else {
        er_neg += wmu_com;
        fr_neg += (wmu_com * mu_com);
      }
    }

    // Loop over angles one last time to set intensities using coefficients
    for (int n=0; n<nang; ++n) {
      Real lab_over_com = LabOverCom(mu[n], beta);

      // Up-going and down-going intensities in the comoving frame
      Real com_ir_pos = (fr_neg * Er_com - er_neg * Fr_com)
                      / (er_pos * fr_neg - fr_pos * er_neg);
      Real com_ir_neg = (er_pos * Fr_com - fr_pos * Er_com)
                      / (er_pos * fr_neg - fr_pos * er_neg);

      // Set lab frame intensities using the transfer coefficient
      Real lorz = sqrt(1.0 / (1.0 - SQR(beta)));
      Real mu_com = (mu[n] - (lorz * beta) * (1.0 - lorz / (lorz + 1.0)
                    * (beta * mu[n]))) / sqrt(lab_over_com);

      if (mu_com > 0.0) {
        lab_ir[n] = com_ir_pos / (SQR(lab_over_com));
      } else {
        lab_ir[n] = com_ir_neg / (SQR(lab_over_com));
      }
    }
  }
} // end anonymous namespace

