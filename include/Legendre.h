// $Header$
//----------------------------------------------------------------------------------------
//                          Legendre
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/14
//
/*
 * @file Legendre.h
 * @brief Computation of Legendre polynomials and their derivatives
 *
 * @details This file contains the implementation of functions for computing Legendre
 * polynomials and their derivatives.
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _LEGENDRE_
#define _LEGENDRE_

#include "matrix.h"
#include <math.h>

/**
 * @brief Computes the Legendre polynomials and their derivatives.
 *
 * @param n Maximum degree of the Legendre polynomials.
 * @param m Maximum order of the Legendre polynomials.
 * @param fi Argument at which the Legendre polynomials are evaluated [rad].
 * @param pnm Output matrix containing the Legendre polynomials.
 * @param dpnm Output matrix containing the derivatives of the Legendre polynomials.
 */
void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm);

#endif