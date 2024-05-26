// $Header$
//----------------------------------------------------------------------------------------
//                          utils
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/23
//
/*
 * @file utils.h
 * @brief Header file for computing the remainder of x divided by y.
 *
 * @param x The dividend.
 * @param y The divisor.
 * @return The remainder of x divided by y.
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _UTILS_
#define _UTILS_

#include <math.h>


/**
 * @brief Computes the remainder of x divided by y.
 *
 * This function computes the remainder of x divided by y.
 * It differs from the standard % operator in handling negative values of x correctly.
 *
 * @param x The dividend.
 * @param y The divisor.
 * @return The remainder of x divided by y.
 */
double custom_mod(double x, double y);

#endif