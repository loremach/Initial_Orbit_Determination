// $Header$
//----------------------------------------------------------------------------------------
//                          G_AccelHarmonic
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/08
//
/*
 * @file G_AccelHarmonic.h
 * @brief Header file for computing the gradient of the Earth's harmonic gravity field
 *
 * @details This header file contains the declaration of the function to compute the gradient of the Earth's harmonic gravity field.
 *
 * @param r Satellite position vector in the true-of-date system
 * @param U Transformation matrix to body-fixed system
 * @param n Gravity model degree
 * @param m Gravity model order
 * @return G Gradient (G = da/dr) in the true-of-date system
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _GACCELHARMONIC_
#define _GACCELHARMONIC_

#include <math.h>
#include "AccelHarmonic.h"

/**
 * @brief Computes the gradient of the Earth's harmonic gravity field.
 *
 * @param r Satellite position vector in the true-of-date system
 * @param U Transformation matrix to body-fixed system
 * @param n Gravity model degree
 * @param m Gravity model order
 * @param G Gradient (G = da/dr) in the true-of-date system
 */
Matrix G_AccelHarmonic(Matrix r, Matrix U, int n_max, int m_max);

#endif