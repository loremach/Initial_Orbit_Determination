// $Header$
//----------------------------------------------------------------------------------------
//                          AccelPointMass
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
* @file AccelPointMass.h
* @brief Header file for computing the perturbational acceleration due to a point mass
 *
 * @details This header file contains the declaration of the function to compute the perturbational acceleration
 * due to a point mass.
 *
 * @param r Satellite position vector
 * @param s Point mass position vector
 * @param GM Gravitational coefficient of the point mass
 * @return a Acceleration (a = d^2r/dt^2)
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _ACCELPOINTMASS_
#define _ACCELPOINTMASS_

#include <math.h>
#include "matrix.h"

/**
 * @brief Computes the perturbational acceleration due to a point mass.
 *
 * @param r Satellite position vector
 * @param s Point mass position vector
 * @param GM Gravitational coefficient of the point mass
 * @param a Acceleration (a = d^2r/dt^2)
 */
Matrix AccelPointMass(Matrix r, Matrix s, double GM);

#endif