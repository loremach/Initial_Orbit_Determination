// $Header$
//----------------------------------------------------------------------------------------
//                          NutAngles
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/28
//
/*
 * @file NutAngles.h
 * @brief Header file for computing the nutation in longitude and obliquity
 *
 * @details This header file contains declarations related to computing the nutation in longitude and obliquity.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @param[out] dpsi Nutation in longitude [rad]
 * @param[out] deps Nutation in obliquity [rad]
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _NUTANGLES_
#define _NUTANGLES_

#include <math.h>
#include "SAT_const.h"

/**
 * @brief Computes the nutation in longitude and obliquity.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @param[out] dpsi Nutation in longitude [rad]
 * @param[out] deps Nutation in obliquity [rad]
 */
void NutAngles (double Mjd_TT, double& dpsi, double& deps);

#endif