// $Source$
//----------------------------------------------------------------------------------------
//                          EqnEquinox
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/02
//
/*
 * @file EqnEquinox.cpp
 * @brief Source file for computing the equation of the equinoxes
 *
 * @details This file contains the declaration of the function to compute the equation of the equinoxes.
 * The equation of the equinoxes dpsi*cos(eps) is the right ascension of the mean equinox referred
 * to the true equator and equinox and is equal to the difference between apparent and mean sidereal time.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return Equation of the equinoxes
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\EqnEquinox.h"

double EqnEquinox(double Mjd_TT)
{
    // Nutation in longitude and obliquity
    double dpsi;
    double deps;
    NutAngles(Mjd_TT, dpsi, deps);

    // Equation of the equinoxes
    return dpsi * cos(MeanObliquity(Mjd_TT));
}