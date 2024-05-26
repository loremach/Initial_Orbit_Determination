// $Header$
//----------------------------------------------------------------------------------------
//                          Frac
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/18
//
/*
 * @file Frac.h
 * @brief Header file for computing the fractional part of a number
 *
 * @details This header file contains the declaration of the function to compute the fractional part of a number.
 *
 * @param x The input number
 * @return The fractional part of the input number (y = x - floor(x))
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _FRAC_
#define _FRAC_

#include "matrix.h"
#include <math.h>

/**
 * @brief Header file for computing the fractional part of a number.
 *
 * @param x The input number
 * @return The fractional part of the input number (y = x - floor(x))
 */
double Frac(double x);

#endif