#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file globals.hpp
//! \brief namespace containing external global variables

namespace Globals {
extern int my_rank, nranks;
}

namespace Uov {

// // python 3.7 --- to generate user output variable enum
//nmu = 16
//enumnames = ["DIVFLX", "ANGFLX", "IMPANG", "CONSTCOEF", "INITFLX", "FINALFLX"]
//text = "GRAV, RADSRC, "
//for m, name in enumerate(enumnames):
//    for n in range(2*nmu):
//        text += f"{name}_{n}, "
//text += "NUM_UOV"
//print(text)

// User output variables
enum {GRAD_P_ACCEL, GRAD_P_ESTIMATE, GRAV_ACCEL, RAD_ACCEL, NET_ACCEL, MDOT, EDOT, CSOUND, TRAD_OVER_TGAS, DIVFLX, SRCFLX, DIV_ER, SRC_ER, DIVFLX_ANG, DIV_ER_ANG, DIVFLX_LAB, SRCFLX_LAB, DIV_ER_LAB, SRC_ER_LAB, DIVFLX_ANG_LAB, DIV_ER_ANG_LAB, TAU_PER_CELL, TAU_FACT, NUM_UOV};
}

#endif // GLOBALS_HPP_
