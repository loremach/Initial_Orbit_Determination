// $Header$
//----------------------------------------------------------------------------------------
//                          EqnEquinox
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/02
//
/*
 * @file EqnEquinox.h
 * @brief Header file for computing the equation of the equinoxes
 *
 * @details This header file contains the declaration of the function to compute the equation of the equinoxes.
 * The equation of the equinoxes dpsi*cos(eps) is the right ascension of the mean equinox referred
 * to the true equator and equinox and is equal to the difference between apparent and mean sidereal time.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return Equation of the equinoxes
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _EQNEQUINOX_
#define _EQNEQUINOX_

#include <math.h>
#include "..\include\NutAngles.h"
#include "..\include\MeanObliquity.h"

/**
 * @brief Computes the equation of the equinoxes.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return Equation of the equinoxes
 */
double EqnEquinox (double Mjd_TT);

#endif