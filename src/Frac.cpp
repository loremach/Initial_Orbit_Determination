// $Source$
//----------------------------------------------------------------------------------------
//                          Frac
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/18
//
/*
 * @file Frac.cpp
 * @brief Computes the fractional part of a number (y = x - floor(x))
 *
 * @details This file contains the implementation of the function to compute the fractional part of a number.
 *
 * @param x The input number
 * @return The fractional part of the input number (y = x - floor(x))
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\Frac.h"

double Frac(double x){
    return x-floor(x);
}