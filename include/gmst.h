// $Header$
//----------------------------------------------------------------------------------------
//                          gmst
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/04
//
/*
 * @file gmst.h
 * @brief Header file for computing the Greenwich Mean Sidereal Time (GMST)
 *
 * @details This header file contains the declaration of the function to compute the Greenwich Mean Sidereal Time.
 *
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return gmstime GMST in [rad]
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _GMST_
#define _GMST_

#include <math.h>
#include "..\include\Frac.h"
#include "..\include\SAT_const.h"

/**
 * @brief Computes the Greenwich Mean Sidereal Time (GMST).
 *
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return gmstime GMST in [rad]
 */
double gmst(double Mjd_UT1);

#endif