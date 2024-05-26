// $Source$
//----------------------------------------------------------------------------------------
//                          gast
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/05
//
/*
 * @file gast.cpp
 * @brief Computes the Greenwich Apparent Sidereal Time (GAST)
 *
 * @details This file contains the implementation of the function to compute the Greenwich Apparent Sidereal Time.
 *
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return gstime GAST in [rad]
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\gast.h"

double gast (double Mjd_UT1){
    return custom_mod (gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), Const::pi2 );
}