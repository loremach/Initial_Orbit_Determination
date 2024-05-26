// $Header$
//----------------------------------------------------------------------------------------
//                          JPL_Eph_DE430
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/10
//
/*
 * @file JPL_Eph_DE430.h
 * @brief Header file for computing the positions of the sun, moon, and nine major planets' equatorial positions using JPL Ephemerides
 *
 * @details This header file contains declarations related to computing the positions of the sun, moon, and nine major planets' equatorial positions using JPL Ephemerides.
 *
 * @param Mjd_TDB Modified Julian Date of TDB
 * @param[out] r_Earth Solar system barycenter (SSB) position
 * @param[out] r_Mars Mars position
 * @param[out] r_Mercury Mercury position
 * @param[out] r_Venus Venus position
 * @param[out] r_Jupiter Jupiter position
 * @param[out] r_Saturn Saturn position
 * @param[out] r_Uranus Uranus position
 * @param[out] r_Neptune Neptune position
 * @param[out] r_Pluto Pluto position
 * @param[out] r_Moon Moon position
 * @param[out] r_Sun Geocentric equatorial position referred to the International Celestial Reference Frame (ICRF)
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _JPL_
#define _JPL_

#include <math.h>
#include "global.h"
#include "Cheb3D.h"

/**
 * @brief Computes the positions of the sun, moon, and nine major planets' equatorial positions using JPL Ephemerides.
 *
 * @param Mjd_TDB Modified Julian Date of TDB
 * @param[out] r_Earth Solar system barycenter (SSB) position
 * @param[out] r_Mars Mars position
 * @param[out] r_Mercury Mercury position
 * @param[out] r_Venus Venus position
 * @param[out] r_Jupiter Jupiter position
 * @param[out] r_Saturn Saturn position
 * @param[out] r_Uranus Uranus position
 * @param[out] r_Neptune Neptune position
 * @param[out] r_Pluto Pluto position
 * @param[out] r_Moon Moon position
 * @param[out] r_Sun Geocentric equatorial position referred to the International Celestial Reference Frame (ICRF)
 */
void JPL_Eph_DE430(Matrix & r_Mercury, Matrix & r_Venus, Matrix & r_Earth, Matrix & r_Mars, 
                    Matrix & r_Jupiter, Matrix & r_Saturn, Matrix & r_Uranus, Matrix & r_Neptune, 
                    Matrix & r_Pluto, Matrix & r_Moon, Matrix & r_Sun, double Mjd_TDB);

#endif