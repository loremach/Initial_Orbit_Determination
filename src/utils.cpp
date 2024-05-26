// $Source$
//----------------------------------------------------------------------------------------
//                          utils
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/23
//
/*
 * @file utils.cpp
 * @brief Computes the remainder of x divided by y.
 *
 * @param x The dividend.
 * @param y The divisor.
 * @return The remainder of x divided by y.
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\utils.h"

double custom_mod(double x, double y){
    return fmod(fmod(x, y) + y, y);
}