// $Header$
//----------------------------------------------------------------------------------------
//                          Accel
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/15
//
/* @file Accel.h
 * @brief Header file for computing the acceleration of an Earth orbiting satellite
 *
 * @details This header file contains the declaration of the function to compute the acceleration of an Earth orbiting satellite due to:
 * - the Earth's harmonic gravity field,
 * - the gravitational perturbations of the Sun and Moon,
 * - the solar radiation pressure, and
 * - the atmospheric drag.
 *
 * @param Mjd_TT Terrestrial Time (Modified Julian Date)
 * @param Y Satellite state vector in the ICRF/EME2000 system
 * @return Acceleration (a = d^2r/dt^2) in the ICRF/EME2000 system
 *
 * @author Lorena Remacha Bordallo
*/
#ifndef _ACCEL_
#define _ACCEL_

#include <math.h>
#include "..\include\Matrix.h"
#include "..\include\IERS.h"
#include "..\include\timediff.h"
#include "..\include\PrecMatrix.h"
#include "..\include\NutMatrix.h"
#include "..\include\PoleMatrix.h"
#include "..\include\GHAMatrix.h"
#include "..\include\Mjday_TDB.h"
#include "..\include\JPL_Eph_DE430.h"
#include "..\include\AccelHarmonic.h"
#include "..\include\AccelPointMass.h"

//----------------------------------------------------------------------------------------
// Matrix Accel(double x, Matrix & Y)
//----------------------------------------------------------------------------------------
/**
 * @brief Computes the acceleration due to various perturbations.
 *
 * @param Mjd_TT Terrestrial Time (Modified Julian Date)
 * @param Y Satellite state vector in the ICRF/EME2000 system
 * @param dY Acceleration (a = d^2r/dt^2) in the ICRF/EME2000 system
*/
//----------------------------------------------------------------------------------------
Matrix Accel(double x, Matrix & Y);

#endif