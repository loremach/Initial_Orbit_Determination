// $Header$
//----------------------------------------------------------------------------------------
//                          AzElPa
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
 * @file AzElPa.h
 * @brief Header file for computing azimuth, elevation, and partial derivatives from local tangent coordinates
 *
 * @details This header file contains the implementation of the function to compute azimuth, elevation,
 * and their partial derivatives with respect to local tangent coordinates in the East-North-Zenith frame.
 *
 * @param s Topocentric local tangent coordinates (East-North-Zenith frame)
 * @param A Azimuth in radians
 * @param E Elevation in radians
 * @param dAds Partials of azimuth with respect to s
 * @param dEds Partials of elevation with respect to s
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _AZELPA_
#define _AZELPA_

#include <math.h>
#include "matrix.h"
#include "SAT_const.h"

/**
 * @brief Computes azimuth, elevation, and partial derivatives from local tangent coordinates.
 *
 * @param s Topocentric local tangent coordinates (East-North-Zenith frame)
 * @param A Azimuth in radians
 * @param E Elevation in radians
 * @param dAds Partials of azimuth with respect to s
 * @param dEds Partials of elevation with respect to s
 */
void AzElPa(Matrix s, double& Az, double& El, Matrix& dAds, Matrix& dEds);

#endif