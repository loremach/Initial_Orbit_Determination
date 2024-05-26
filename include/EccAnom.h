// $Header$
//----------------------------------------------------------------------------------------
//                          EccAnom
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
 * @file EccAnom.h
 * @brief Header file for computing the eccentric anomaly for elliptic orbits
 *
 * @details This header file contains the declaration of the function to compute the eccentric anomaly
 * for elliptic orbits given the mean anomaly and the eccentricity.
 *
 * @param M Mean anomaly in radians
 * @param e Eccentricity of the orbit (in the range [0,1])
 * @return Eccentric anomaly in radians
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _ECCANOM_
#define _ECCANOM_

#include <math.h>
#include "SAT_const.h"

/**
 * @brief Computes the eccentric anomaly for elliptic orbits.
 *
 * @param M Mean anomaly in radians
 * @param e Eccentricity of the orbit (in the range [0,1])
 * @return Eccentric anomaly in radians
 */
double EccAnom (double M, double e);

#endif