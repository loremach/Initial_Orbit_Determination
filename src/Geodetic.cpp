// $Source$
//----------------------------------------------------------------------------------------
//                          Geodetic
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
 * @file Geodetic.cpp
 * @brief Computes geodetic coordinates (Longitude [rad], latitude [rad], altitude [m]) from a given position vector (r [m])
 *
 * @details This file contains the implementation of the function to compute geodetic coordinates from a given position vector.
 *
 * @param r Position vector in meters
 * @param[out] lon Longitude in radians
 * @param[out] lat Latitude in radians
 * @param[out] alt Altitude in meters
 *
 * @author Lorena Remacha Bordallo
*/
#include "..\include\Geodetic.h"

void Geodetic(Matrix r, double& lon, double& lat, double& h){
    double R_equ = Const::R_Earth;
    double f     = Const::f_Earth;

    double epsRequ = 2.22044604925031e-16*R_equ;        // Convergence criterion
    double e2      = f*(2.0-f);        // Square of eccentricity

    double X = r(1);                   // Cartesian coordinates
    double Y = r(2);
    double Z = r(3);
    double rho2 = X*X + Y*Y;           // Square of distance from z-axis

    // Check validity of input data
    if (norm(r)==0.0){
        cout << "invalid input in Geodetic constructor\n";
        exit(EXIT_FAILURE);
        lon = 0.0;
        lat = 0.0;
        h   = - R_equ;
    }

    // Iteration 
    double dZ = e2*Z;
    double ZdZ, Nh, SinPhi, N, dZ_new;
    while(true){
        ZdZ    =  Z + dZ;
        Nh     =  sqrt ( rho2 + ZdZ*ZdZ ); 
        SinPhi =  ZdZ / Nh;                    // Sine of geodetic latitude
        N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
        dZ_new =  N*e2*SinPhi;

        if (fabs(dZ-dZ_new) < epsRequ )
            break;
        
        dZ = dZ_new;
    }

    // Longitude, latitude, altitude
    lon = atan2 ( Y, X );
    lat = atan2 ( ZdZ, sqrt(rho2) );
    h   = Nh - N;
}