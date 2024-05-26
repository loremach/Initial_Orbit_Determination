// $Source$
//----------------------------------------------------------------------------------------
//                          sign
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/18
//
/*
 * @file sign.cpp
 * @brief Returns the absolute value of a with the sign of b
 *
 * @details This file contains the implementation of the function to return the absolute value of a with the sign of b.
 *
 * @param a Input value
 * @param b Reference value
 * @return Absolute value of a with the sign of b
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\sign_.h"

double sign_(double a, double b){
    double result;
    if (b>=0.0)
        result = fabs(a);
    else
        result = - fabs(a);

    return result;
}