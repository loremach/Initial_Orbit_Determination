// $Header$
//----------------------------------------------------------------------------------------
//                          elements
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/05
//
/*
 * @file elements.h
 * @brief Header file for computing the osculating Keplerian elements from the satellite state vector for elliptic orbits
 *
 * @details This header file contains the declaration of the function to compute the osculating Keplerian elements
 * from the satellite state vector for elliptic orbits. The function cannot be used with state vectors
 * describing a circular or non-inclined orbit.
 *
 * @param y State vector (x, y, z, vx, vy, vz)
 * @param p Semilatus rectum in meters
 * @param a Semimajor axis
 * @param e Eccentricity
 * @param i Inclination in radians
 * @param Omega Longitude of the ascending node in radians
 * @param omega Argument of pericenter in radians
 * @param M Mean anomaly in radians
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _ELEMENTS_
#define _ELEMENTS_

#include <math.h>
#include "matrix.h"
#include "SAT_const.h"

/**
 * @brief Computes the osculating Keplerian elements from the satellite state vector for elliptic orbits.
 *
 * @param y State vector (x, y, z, vx, vy, vz)
 * @param p Semilatus rectum in meters
 * @param a Semimajor axis
 * @param e Eccentricity
 * @param i Inclination in radians
 * @param Omega Longitude of the ascending node in radians
 * @param omega Argument of pericenter in radians
 * @param M Mean anomaly in radians
 */
void elements(double & p, double & a, double & e, double & i, double & Omega, double & omega, double & M, Matrix y);

#endif