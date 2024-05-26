// $Header$
//----------------------------------------------------------------------------------------
//                          Position
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
 * @file Position.h
 * @brief Header file for computing the position vector (r [m]) from geodetic coordinates
 *
 * @details This header file contains declarations related to computing the position vector (r [m]) from geodetic coordinates.
 *
 * @param lon Geodetic East longitude [rad]
 * @param lat Geodetic latitude [rad]
 * @param alt Altitude [m]
 * @param[out] r Position vector (r [m])
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _POSITION_
#define _POSITION_

#include "matrix.h"
#include <math.h>
#include "SAT_Const.h"

/**
 * @brief Computes the position vector (r [m]) from geodetic coordinates.
 *
 * @param lon Geodetic East longitude [rad]
 * @param lat Geodetic latitude [rad]
 * @param alt Altitude [m]
 * @param[out] r Position vector (r [m])
 */
Matrix Position(double lon, double lat, double h);

#endif