// $Header$
//----------------------------------------------------------------------------------------
//                          Geodetic
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
 * @file Geodetic.h
 * @brief Header file for computing geodetic coordinates (Longitude [rad], latitude [rad], altitude [m]) from a given position vector (r [m])
 *
 * @details This header file contains the declaration of the function to compute geodetic coordinates from a given position vector.
 *
 * @param r Position vector in meters
 * @param[out] lon Longitude in radians
 * @param[out] lat Latitude in radians
 * @param[out] alt Altitude in meters
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _GEODETIC_
#define _GEODETIC_

#include <math.h>
#include "SAT_Const.h"
#include "matrix.h"

/**
 * @brief Computes geodetic coordinates (Longitude [rad], latitude [rad], altitude [m]) from a given position vector (r [m]).
 *
 * @param r Position vector in meters
 * @param[out] lon Longitude in radians
 * @param[out] lat Latitude in radians
 * @param[out] alt Altitude in meters
 */
void Geodetic(Matrix r, double& lon, double& lat, double& h);

#endif