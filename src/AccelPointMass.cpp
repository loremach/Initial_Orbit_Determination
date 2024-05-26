// $Source$
//----------------------------------------------------------------------------------------
//                          AccelPointMass
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
* @file AccelPointMass.cpp
* @brief Computes the perturbational acceleration due to a point mass
 *
 * @details This file contains the declaration of the function to compute the perturbational acceleration
 * due to a point mass.
 *
 * @param r Satellite position vector
 * @param s Point mass position vector
 * @param GM Gravitational coefficient of the point mass
 * @return a Acceleration (a = d^2r/dt^2)
 *
 * @author Lorena Remacha Bordallo
*/
#include "..\include\AccelPointMass.h"


Matrix AccelPointMass(Matrix r, Matrix s, double GM){
    // Relative position vector of satellite w.r.t. point mass 
    Matrix d = r - s;
    // Acceleration 
    return (d/pow(norm(d),3) + s/pow(norm(s),3)) * (-GM);
}