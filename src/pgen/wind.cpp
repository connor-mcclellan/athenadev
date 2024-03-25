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

// file-scope variables for MESA profile data
AthenaArray<Real> mesa_og;
AthenaArray<Real> mesa_in;

void SetLabIr(int nang, Real beta, Real Er, Real Fr, int k, int j, int i,
              Real *lab_ir, Real *mu, Real *wmu);

bool is_between(Real x, Real a, Real b);
bool is_between(Real x, Real a, Real b, Real &t);

// MESA variable data fields
enum {RHO, GRAV, RADIUS, VELOCITY, PRAD, PGAS, MU, OPACITY, LUMINOSITY, TAU,
      NUM_MESA};

// Dimensional quantities
Real r0;        // Length
Real v0;        // Velocity
Real rho0;      // Density
Real temp0;     // Temperature
Real temp04;
Real kappa0;    // Opacity
Real mmw0;

// Physical constants
Real clight;
Real arad;
Real kb;
Real mp;

Real energy_alpha;

} // namespace

void RadInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt, int is, int ie, int js, int je, int ks,
                int ke, int ngh);

void HydroInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt, int il, int iu, int jl,
                  int ju, int kl, int ku, int ngh);

void UserOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);


void Mesh::InitUserMeshData(ParameterInput *pin) {
  Real crat = pin->GetReal("radiation", "crat");
  Real prat = pin->GetReal("radiation", "prat");

  Real GM = pin->GetReal("problem", "GM");
  Real mass = pin->GetReal("problem", "mass");
  Real GM_cgs = mass * Constants::solar_mass_cgs * Constants::grav_const_cgs;
  energy_alpha = pin->GetReal("problem", "energy_alpha");

  clight = Constants::speed_of_light_cgs;
  kb = Constants::k_boltzmann_cgs;
  arad = Constants::radiation_aconst_cgs;
  mp = Constants::hydrogen_mass_cgs;

  v0 = clight / crat;
  mmw0 = 4./3;
  temp0 = SQR(v0) * mmw0 * mp / kb;
  temp04 = std::pow(temp0, 4.0);
  r0 = GM_cgs/SQR(clight);
  kappa0 = 0.2;
  rho0 = 1. / kappa0 / r0;

  EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, RadInnerX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, HydroInnerX1);
  return;
} // InitUserMeshData


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {

  pnrrad->EnrollOpacityFunction(UserOpacity);

  // ATHENA USER OUTPUT I/O
  // ======================

  // Get the size of the mesh, including ghost zones
  int nx1 = pmy_mesh->mesh_size.nx1 + 2 * NGHOST;
  int nx2 = pmy_mesh->mesh_size.nx2 + 2 * NGHOST;
  int nx3 = pmy_mesh->mesh_size.nx3 + 2 * NGHOST;

  AllocateRealUserMeshBlockDataField(Uov::NUM_UOV);
  for (int l=0; l<Uov::NUM_UOV; l++) {
    ruser_meshblock_data[l].NewAthenaArray(nx3, nx2, nx1);
  }
  AllocateUserOutputVariables(Uov::NUM_UOV);

  // Output variables that are not angle dependent - one output per timestep
  SetUserOutputVariableName(Uov::CSOUND, "csound");
  SetUserOutputVariableName(Uov::GRAD_P, "grad_P");
  SetUserOutputVariableName(Uov::ANGFLX_SUM, "angflx_sum");
  SetUserOutputVariableName(Uov::DIVFLX_SUM, "divflx_sum");
  SetUserOutputVariableName(Uov::INITFLX_SUM, "initflx_sum");
  SetUserOutputVariableName(Uov::TAU, "tau");
  SetUserOutputVariableName(Uov::GRAVSRC, "gravsrc");
  SetUserOutputVariableName(Uov::RADSRC, "radsrc");

  // MESA DATA I/O
  // =============

  // MESA profile data temporary arrays
  std::vector<std::string> headers, hdata, vars;
  std::vector<Real> vdata;

  // Load data from MESA profile into temporary arrays
  std::string filename = pin->GetString("problem", "profile");
  MesaReader(filename.c_str(), headers, hdata, vars, vdata);
  int nmesaradii = vdata.size() / vars.size(); // (nrows * nvars) / (nvars)

  // Get header information from MESA for unit conversion
  Real lsun, msun;
  GetMesaHeader("lsun", headers, hdata, lsun);
  GetMesaHeader("msun", headers, hdata, msun);

  // Copy variable data from temporary arrays to AthenaArray
  mesa_og.NewAthenaArray(NUM_MESA, nmesaradii);
  GetMesaData("logRho", vars, vdata, &mesa_og(RHO, 0));            // cell center
  GetMesaData("grav", vars, vdata, &mesa_og(GRAV, 0));             // outer face (?)
  GetMesaData("radius_cm", vars, vdata, &mesa_og(RADIUS, 0));      // outer face
  GetMesaData("velocity", vars, vdata, &mesa_og(VELOCITY, 0));     // outer face
  GetMesaData("prad", vars, vdata, &mesa_og(PRAD, 0));             // cell center
  GetMesaData("pgas", vars, vdata, &mesa_og(PGAS, 0));             // cell center
  GetMesaData("mu", vars, vdata, &mesa_og(MU, 0));                 // cell center (?)
  GetMesaData("log_opacity", vars, vdata, &mesa_og(OPACITY, 0));   // cell center
  GetMesaData("logL", vars, vdata, &mesa_og(LUMINOSITY, 0));       // outer face
  GetMesaData("tau", vars, vdata, &mesa_og(TAU, 0));               // cell center

  // Convert MESA variables to code units
  for (int l=0; l<NUM_MESA; ++l) {
    for (int n=0; n<nmesaradii; n++) {
      if (l==RHO)
        mesa_og(l,n) = std::pow(10.0, mesa_og(l,n)) / rho0;
      else if (l==GRAV)
        mesa_og(l,n) /= r0/SQR(v0);
      else if (l==RADIUS)
        mesa_og(l,n) /= r0;
      else if (l==VELOCITY)
        mesa_og(l,n) /= v0;
      else if (l==PRAD)
        mesa_og(l,n) /= arad * temp04;
      else if (l==PGAS)
        mesa_og(l,n) /= rho0*SQR(v0);
      else if (l==OPACITY)
        mesa_og(l,n) = std::pow(10.0, mesa_og(l,n)) / kappa0;
      else if (l==LUMINOSITY)
        mesa_og(l,n) = std::pow(10.0, mesa_og(l,n)) * lsun /
                       (4.0 * PI * SQR(r0) * clight * arad * temp04);
    }
  }

  // Interpolate the mesa variables onto the Athena++ cell-centered grid
  mesa_in.NewAthenaArray(NUM_MESA, nx1);
  for (int i=is-NGHOST; i<=ie+NGHOST; ++i) {
    int n;
    Real fraction = 0.0;
    Real r = pcoord->x1v(i);
    if (!is_between(r, mesa_og(RADIUS,nmesaradii-1), mesa_og(RADIUS,0))) {
      std::stringstream msg;
      msg << "### FATAL ERROR in MeshBlock::InitUserMeshBlockData" << std::endl
          << "Attempted MESA interpolation out of bounds (r=" << r << " not in ["
          << mesa_og(RADIUS,nmesaradii-1) << ", " << mesa_og(RADIUS,0) << "])."
          << std::endl;
      ATHENA_ERROR(msg);
    } else {
      for (n=nmesaradii-1; n>0; --n) {
        bool bracketed = is_between(r, mesa_og(RADIUS,n), mesa_og(RADIUS,n-1), fraction);
        if (bracketed) {
          break;
        }
      }
    }
    for (int l=0; l<NUM_MESA; ++l) {
      mesa_in(l,i) = (fraction)*mesa_og(l,n-1) + (1.0-fraction)*mesa_og(l,n);
    }
  }
  mesa_og.DeleteAthenaArray();
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
              user_out_var(Uov::GRAD_P,k,j,i) = diff_P / pcoord->dx1f(i); // Pressure gradient force
            } else {
              user_out_var(Uov::GRAD_P,k,j,i) = 0.0;
            }
          } else if (l==Uov::CSOUND) {
            Real rho = phydro->u(IDN, k, j, i);
            Real press = phydro->w(IPR,k,j,i);
            user_out_var(Uov::CSOUND,k,j,i) = std::sqrt(gamma * press / rho);
          } else if (l==Uov::GRAVSRC) {
            Real rho = phydro->u(IDN, k, j, i);
            Real mx = phydro->u(IM1, k, j, i);
            Real my = phydro->u(IM2, k, j, i);
            Real mz = phydro->u(IM3, k, j, i);

            // Left and right cell faces
            Real rl = pcoord->x1f(i);
            Real rr = pcoord->x1f(i + 1);
            Real rc = pcoord->x1v(i);

            // Gravitational parameter
            Real gm = SQR(pnrrad->crat);

            // Gravitational potential
            Real phil = -gm/rl;
            Real phir = -gm/rr;

            // Source term - force per unit volume
            Real dt = pmy_mesh->dt;
            Real src = - dt * rho * (phir - phil) / (rr - rl);

            // Momentum change per timestep due to gravity
            user_out_var(Uov::GRAVSRC, k, j, i) = src / dt;

          } else if ((l==Uov::TAU) && (i != ie+NGHOST)) {
            Real tau_avg = 0.5*(pnrrad->sigma_s(k,j,i,0) + pnrrad->sigma_s(k,j,i+1,0)) * pcoord->dx1v(i);
            user_out_var(l,k,j,i) = tau_avg;
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

  // Initialize the gas and radiation field
  for (int k = ks; k <= ke; k++) {
    for (int j = js; j <= je; j++) {
      for (int i = is - NGHOST; i <= ie + NGHOST; i++) {

        Real x1v = pcoord->x1v(i);
        Real rho = mesa_in(RHO, i);
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = mesa_in(VELOCITY, i) * rho;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        Real Fr = mesa_in(LUMINOSITY, i) / SQR(x1v);
        //Real Er = 3.0 * mesa_in(PRAD, i);
        //Real temp = std::pow(Er, 0.25);
        Real temp = mesa_in(PGAS, i) / rho;
        Real Er = std::pow(temp, 4.0);

        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = mesa_in(PGAS, i) / (gamma - 1.0);
          phydro->u(IEN, k, j, i) += 0.5 * SQR(phydro->u(IM1, k, j, i)) / rho;
          phydro->u(IEN, k, j, i) += 0.5 * SQR(phydro->u(IM2, k, j, i)) / rho;
          phydro->u(IEN, k, j, i) += 0.5 * SQR(phydro->u(IM3, k, j, i)) / rho;
        }

        if (IM_RADIATION_ENABLED || NR_RADIATION_ENABLED) {
            Real *lab_ir = &(pnrrad->ir(k, j, i, 0));
            Real *mu     = &(pnrrad->mu(0, k, j, i, 0));
            Real *wmu    = &(pnrrad->wmu(0));
            Real beta = phydro->u(IM1,k,j,i)/rho/pnrrad->crat;
            SetLabIr(pnrrad->nang, beta, Er, Fr, k, j, i, lab_ir, mu, wmu);
        }
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

          Real kappa_ff = 2.87e24 * (rho * rho0) / pow(temp * temp0, 3.5);
          kappa_ff /= kappa0;

          // Scattering opacity
          //Real kappa_s = 1. / (1. + std::pow(temp * temp0 / 4.5e8, 0.86));
          Real kappa_s = 1.0;

          // Set all the matter-radiation coupling coefficients
          pnrrad->sigma_s(k, j, i, ifr) = rho * kappa_s;
          pnrrad->sigma_a(k, j, i, ifr) = rho * kappa_ff * 0.1;
          pnrrad->sigma_p(k, j, i, ifr) = rho * kappa_ff;
          pnrrad->sigma_pe(k, j, i, ifr) = rho * kappa_ff;
        }
      }
    }
  }

}

void HydroInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt, int il, int iu, int jl,
                  int ju, int kl, int ku, int ngh) {

  NRRadiation *pnrrad = pmb->pnrrad;
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {

      // Calculate energy density in first active zone
      Real Er_above = 0.0;
      for (int n=0; n < pnrrad->nang; n++) {
        Er_above += pnrrad->ir(k,j,il,n) * pnrrad->wmu(n);
      }

      for (int i=il-1; i>=il-ngh; --i) {
        Real x1v = pco->x1v(i);
        Real Fr = mesa_in(LUMINOSITY, i) / SQR(x1v);
        Real tau_avg = 0.5*(pnrrad->sigma_s(k,j,i,0) + pnrrad->sigma_s(k,j,i+1,0)) * pco->dx1v(i);
        Real Er = Er_above + energy_alpha * 3.0 * Fr * tau_avg;
        Real temp = std::pow(Er, 0.25);

        for (int n=0; n<=NHYDRO; n++) {
          if (n == IDN) {
            prim(n, k, j, i) = mesa_in(RHO, i);
          } else if (n == IPR) {
            prim(n, k, j, i) = mesa_in(RHO, i) * temp;
          } else if (n == IVX) {
            prim(n, k, j, i) = mesa_in(VELOCITY, i);
          } else {
            prim(n, k, j, i) = prim(n, k, j, il);
          }
        }
        Er_above = Er;
      }
    }
  }
  return;
}

void RadInnerX1(MeshBlock *pmb, Coordinates *pco, NRRadiation *pnrrad,
                const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &ir,
                Real time, Real dt, int is, int ie, int js, int je, int ks,
                int ke, int ngh) {

  Real gamma = pmb->peos->GetGamma();

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {

      // Calculate energy density in first active zone
      Real Er_above = 0.0;
      for (int n=0; n < pnrrad->nang; n++) {
        Er_above += ir(k,j,is,n) * pnrrad->wmu(n);
      }

      for (int i=is-1; i>=is-ngh; --i) {
        Real x1v = pco->x1v(i);
        Real Fr = mesa_in(LUMINOSITY, i) / SQR(x1v);
        Real tau_avg = 0.5*(pnrrad->sigma_s(k,j,i,0) + pnrrad->sigma_s(k,j,i+1,0)) * pco->dx1v(i);
        Real Er = Er_above + energy_alpha * 3.0 * Fr * tau_avg;

        Real *lab_ir = &(ir(k, j, i, 0));
        Real *mu     = &(pnrrad->mu(0, k, j, i, 0));
        Real *wmu    = &(pnrrad->wmu(0));
        Real beta = w(IVX,k,j,i)/pnrrad->crat;
        SetLabIr(pnrrad->nang, beta, Er, Fr, k, j, i, lab_ir, mu, wmu);

        Er_above = Er;
      }
    }
  }
}



namespace {

bool is_between(Real x, Real a, Real b) {
  if (std::signbit(a - x) != std::signbit(b - x)) {
    return true;
  } else {
    return false;
  }
}

bool is_between(Real x, Real a, Real b, Real &t) {
  if (std::signbit(a - x) != std::signbit(b - x)) {
    t = (x - a) / (b - a);
    return true;
  } else {
    return false;
  }
}

void SetLabIr(int nang, Real beta, Real Er, Real Fr, int k, int j, int i,
              Real *lab_ir, Real *mu, Real *wmu) {

  Real er_pos, fr_pos;
  Real er_neg, fr_neg;

  // Transformation coefficients
  // tran_coef is (\nu_0/\nu) in Mihalas - comoving frequency over lab
  // frequency
  Real lorz = std::sqrt(1.0 / (1.0 - SQR(beta)));
  Real tran_coef[nang];
  Real wmu_cm[nang];
  Real mu_cm[nang];
  Real cm_to_lab[nang];

  // Transform weights and angles to comoving frame
  Real wnorm_com = 0.;
  for (int n = 0; n < nang; ++n) {
    Real vnc = 1.0 - beta * mu[n];
    tran_coef[n] = lorz * vnc;
    cm_to_lab[n] = tran_coef[n] * tran_coef[n] * tran_coef[n] * tran_coef[n];
    mu_cm[n] =
        (1.0 / tran_coef[n]) *
        (mu[n] - lorz * beta * (1.0 - lorz * mu[n] * beta / (lorz + 1.0)));
    wmu_cm[n] = wmu[n] / (tran_coef[n] * tran_coef[n]);
    wnorm_com += wmu_cm[n];
  }

  // Ensure weights sum to 1 in comoving frame
  for (int n = 0; n < nang; ++n) {
    wmu_cm[n] /= wnorm_com;
  }

  // Sum upward and downward intensity coefficients
  er_pos = 0., fr_pos = 0.;
  er_neg = 0., fr_neg = 0.;
  for (int n = 0; n < nang; ++n) {
    if (mu_cm[n] > 0.0) {
      er_pos += wmu_cm[n];
      fr_pos += (wmu_cm[n] * mu_cm[n]);
    } else {
      er_neg += wmu_cm[n];
      fr_neg += (wmu_cm[n] * mu_cm[n]);
    }
  }

  Real ir_pos = (fr_neg * Er - er_neg * Fr) / (er_pos * fr_neg - fr_pos * er_neg);
  Real ir_neg = (er_pos * Fr - fr_pos * Er) / (er_pos * fr_neg - fr_pos * er_neg);

  //printf("zone: %d   ir_pos: %g   ir_neg: %g   cm_to_lab: %g\n", i, ir_pos,
  //       ir_neg, cm_to_lab[0]);

  // set intensities using coefficients
  for (int n = 0; n < nang; ++n) {
    if (mu_cm[n] >= 0.0) {
      lab_ir[n] = ir_pos / cm_to_lab[n];
    } else {
      lab_ir[n] = ir_neg / cm_to_lab[n];
    }
  }
}
} // namespace
