// $Header$
//----------------------------------------------------------------------------------------
//                          Cheb3D
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
 * @file Cheb3D.h
 * @brief Header file for Chebyshev approximation of 3-dimensional vectors
 *
 * @details This header file contains the declaration of the function for Chebyshev approximation
 * of 3-dimensional vectors.
 *
 * @param N Number of coefficients
 * @param Ta Begin interval
 * @param Tb End interval
 * @param Cx Coefficients of Chebyshev polynomial (x-coordinate)
 * @param Cy Coefficients of Chebyshev polynomial (y-coordinate)
 * @param Cz Coefficients of Chebyshev polynomial (z-coordinate)
 * @param t The evaluation point within the interval [Ta, Tb]
 * @param result The resulting 3-dimensional vector approximation
 *
 * @author Lorena Remacha Bordallo
*/
#ifndef _CHEB3D_
#define _CHEB3D_

#include <math.h>
#include "matrix.h"

/**
 * @brief Computes the Chebyshev approximation of 3-dimensional vectors.
 *
 * @param N Number of coefficients
 * @param Ta Begin interval
 * @param Tb End interval
 * @param Cx Coefficients of Chebyshev polynomial (x-coordinate)
 * @param Cy Coefficients of Chebyshev polynomial (y-coordinate)
 * @param Cz Coefficients of Chebyshev polynomial (z-coordinate)
 * @param t The evaluation point within the interval [Ta, Tb]
 * @param result The resulting 3-dimensional vector approximation
 */
Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix Cx, Matrix Cy, Matrix Cz);

#endif