// $Source$
//----------------------------------------------------------------------------------------
//                          AzElPa
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
 * @file AzElPa.cpp
 * @brief Computes azimuth, elevation, and partial derivatives from local tangent coordinates
 *
 * @details This file contains the implementation of the function to compute azimuth, elevation,
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


#include "..\include\AzElPa.h"

void AzElPa(Matrix s, double& Az, double& El, Matrix& dAds, Matrix& dEds){

    double rho = sqrt(s(1)*s(1)+s(2)*s(2));

    // Angles
    Az = atan2(s(1),s(2));

    if (Az<0.0) {
        Az = Az+Const::pi2;
    }

    El = atan ( s(3) / rho );

    // Partials
    dAds(1) = s(2)/(rho*rho);    dAds(2) = -s(1)/(rho*rho);     dAds(3) = 0.0;
    dEds(1) = -s(1)*s(3)/rho;    dEds(2) = -s(2)*s(3)/rho;      dEds(3) = rho;
    dEds = dEds / dot(s,s);
}