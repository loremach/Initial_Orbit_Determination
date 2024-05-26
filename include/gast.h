// $Header$
//----------------------------------------------------------------------------------------
//                          gast
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/05
//
/*
 * @file gast.h
 * @brief Header file for computing the Greenwich Apparent Sidereal Time (GAST)
 *
 * @details This header file contains the declaration of the function to compute the Greenwich Apparent Sidereal Time.
 *
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return gstime GAST in [rad]
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _GAST_
#define _GAST_

#include <math.h>
#include "EqnEquinox.h"
#include "gmst.h"
#include "utils.h"
#include "SAT_const.h"

/**
 * @brief Computes the Greenwich Apparent Sidereal Time (GAST).
 *
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return gstime GAST in [rad]
 */
double gast (double Mjd_UT1);

#endif