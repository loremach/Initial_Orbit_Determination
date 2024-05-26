// $Source$
//----------------------------------------------------------------------------------------
//                          gmst
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/04
//
/*
 * @file gmst.cpp
 * @brief Computes the Greenwich Mean Sidereal Time (GMST)
 *
 * @details This file contains the implementation of the function to compute the Greenwich Mean Sidereal Time.
 *
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return gmstime GMST in [rad]
 *
 * @author Lorena Remacha Bordallo
*/
#include "..\include\gmst.h"

double gmst(double Mjd_UT1){
    double Secs = 86400.0;                       // Seconds per day
    double MJD_J2000 = 51544.5;

    double Mjd_0 = floor(Mjd_UT1);
    double UT1   = Secs*(Mjd_UT1-Mjd_0);         // [s]
    double T_0   = (Mjd_0  -MJD_J2000)/36525.0;
    double T     = (Mjd_UT1-MJD_J2000)/36525.0;

    double gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1 + (0.093104-6.2e-6*T)*T*T;    // [s]

    return 2*Const::pi*Frac(gmst/Secs);       // [rad], 0..2pi
}