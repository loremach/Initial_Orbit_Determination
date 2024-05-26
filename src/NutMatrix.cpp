// $Source$
//----------------------------------------------------------------------------------------
//                          NutMatrix
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/02
//
/*
 * @file NutMatrix.cpp
 * @brief Transformation from mean to true equator and equinox
 *
 * @details This file contains the implementation of the function to perform the transformation from mean to true equator and equinox.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @param[out] NutMat Nutation matrix
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\NutMatrix.h"

Matrix NutMatrix(double Mjd_TT)
{
    // Mean obliquity of the ecliptic
    double eps = MeanObliquity(Mjd_TT);

    double dpsi;
    double deps;
    // Nutation in longitude and obliquity
    NutAngles(Mjd_TT, dpsi, deps);

    // Transformation from mean to true equator and equinox
 
    Matrix r1 = R_x(-eps - deps) ;
    Matrix r2 =  R_z(-dpsi);
    Matrix r3 = R_x(+eps);

    return  r1*r2*r3;

}