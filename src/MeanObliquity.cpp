// $Source$
//----------------------------------------------------------------------------------------
//                          MeanObliquity
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
 * @file MeanObliquity.cpp
 * @brief Computes the mean obliquity of the ecliptic
 *
 * @details This file contains the implementation of the function to compute the mean obliquity of the ecliptic.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return MOblq Mean obliquity of the ecliptic [rad]
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\MeanObliquity.h"

double MeanObliquity (double Mjd_TT){
    double T = (Mjd_TT-Const::MJD_J2000)/36525.0;

    return Const::Rad *( 84381.448/3600.0-(46.8150+(0.00059-0.001813*T)*T)*T/3600.0);
}