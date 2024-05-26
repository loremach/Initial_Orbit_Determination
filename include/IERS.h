// $Header$
//----------------------------------------------------------------------------------------
//                          IERS
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/24
//
/*
 * @file IERS.h
 * @brief Header file for management of IERS time and polar motion data
 *
 * @details This header file contains declarations related to the management of IERS time and polar motion data.
 *
 * @param x_pole Output parameter for x-coordinate of the Celestial Intermediate Pole (CIP) [rad].
 * @param y_pole Output parameter for y-coordinate of the Celestial Intermediate Pole (CIP) [rad].
 * @param UT1_UTC Output parameter for the difference between Universal Time (UT1) and Coordinated Universal Time (UTC) [s].
 * @param LOD Output parameter for the length of day [s].
 * @param dpsi Output parameter for the nutation correction in longitude [rad].
 * @param deps Output parameter for the nutation correction in obliquity [rad].
 * @param dx_pole Output parameter for the rate of change of the x-coordinate of the Celestial Intermediate Pole (CIP) [rad/s].
 * @param dy_pole Output parameter for the rate of change of the y-coordinate of the Celestial Intermediate Pole (CIP) [rad/s].
 * @param TAI_UTC Output parameter for the difference between International Atomic Time (TAI) and Coordinated Universal Time (UTC) [s].
 * @param eop Matrix containing Earth Orientation Parameters (EOP) data.
 * @param Mjd_UTC Modified Julian Date (UTC).
 * @param interp Interpolation type ('l' for linear, 'n' for nearest neighbor).
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _IERS_
#define _IERS_

#include "matrix.h"
#include <math.h>
#include "SAT_Const.h"

/**
 * @brief Function to manage IERS time and polar motion data.
 *
 * @param x_pole Output parameter for x-coordinate of the Celestial Intermediate Pole (CIP) [rad].
 * @param y_pole Output parameter for y-coordinate of the Celestial Intermediate Pole (CIP) [rad].
 * @param UT1_UTC Output parameter for the difference between Universal Time (UT1) and Coordinated Universal Time (UTC) [s].
 * @param LOD Output parameter for the length of day [s].
 * @param dpsi Output parameter for the nutation correction in longitude [rad].
 * @param deps Output parameter for the nutation correction in obliquity [rad].
 * @param dx_pole Output parameter for the rate of change of the x-coordinate of the Celestial Intermediate Pole (CIP) [rad/s].
 * @param dy_pole Output parameter for the rate of change of the y-coordinate of the Celestial Intermediate Pole (CIP) [rad/s].
 * @param TAI_UTC Output parameter for the difference between International Atomic Time (TAI) and Coordinated Universal Time (UTC) [s].
 * @param eop Matrix containing Earth Orientation Parameters (EOP) data.
 * @param Mjd_UTC Modified Julian Date (UTC).
 * @param interp Interpolation type ('l' for linear, 'n' for nearest neighbor).
 */
void IERS(double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi,
          double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC, Matrix eop, 
          double Mjd_UTC, char interp = 'n');


#endif