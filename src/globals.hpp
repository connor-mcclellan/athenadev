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
enum {TAU, GRAV, RADSRC, DIVFLX_0, DIVFLX_1, DIVFLX_2, DIVFLX_3, DIVFLX_4, DIVFLX_5, DIVFLX_6, DIVFLX_7, ANGFLX_0, ANGFLX_1, ANGFLX_2, ANGFLX_3, ANGFLX_4, ANGFLX_5, ANGFLX_6, ANGFLX_7, IMPANG_0, IMPANG_1, IMPANG_2, IMPANG_3, IMPANG_4, IMPANG_5, IMPANG_6, IMPANG_7, CONSTCOEF_0, CONSTCOEF_1, CONSTCOEF_2, CONSTCOEF_3, CONSTCOEF_4, CONSTCOEF_5, CONSTCOEF_6, CONSTCOEF_7, INITFLX_0, INITFLX_1, INITFLX_2, INITFLX_3, INITFLX_4, INITFLX_5, INITFLX_6, INITFLX_7, FINALFLX_0, FINALFLX_1, FINALFLX_2, FINALFLX_3, FINALFLX_4, FINALFLX_5, FINALFLX_6, FINALFLX_7, NUM_UOV};
}

#endif // GLOBALS_HPP_
