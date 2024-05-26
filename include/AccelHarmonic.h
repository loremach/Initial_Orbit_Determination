// $Header$
//----------------------------------------------------------------------------------------
//                          AccelHarmonic
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/05
//
/*
* @file AccelHarmonic.h
* @brief Header file for computing the acceleration due to the harmonic gravity field of the central body
*
* @details This header file contains the implementation of the function to compute the acceleration
        * due to the harmonic gravity field of the central body.
*
* @param r Satellite position vector in the inertial system
* @param E Transformation matrix to body-fixed system
* @param n_max Maximum degree
* @param m_max Maximum order (m_max <= n_max; m_max = 0 for zonals only)
* @return a Acceleration (a = d^2r/dt^2)
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _ACCELHARMONIC_
#define _ACCELHARMONIC_

#include <math.h>
#include "matrix.h"
#include "Legendre.h"
#include "global.h"

/**
 * @brief Computes the acceleration due to the harmonic gravity field of the central body.
 *
 * @param r Satellite position vector in the inertial system
 * @param E Transformation matrix to body-fixed system
 * @param n_max Maximum degree
 * @param m_max Maximum order (m_max <= n_max; m_max = 0 for zonals only)
 * @param a Acceleration (a = d^2r/dt^2)
 */
Matrix AccelHarmonic(Matrix r, Matrix E, int n_max, int m_max);

#endif