// $Source$
//----------------------------------------------------------------------------------------
//                          PoleMatrix
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/03
//
/*
 * @file PoleMatrix.cpp
 * @brief Transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date
 *
 * @details This file contains the implementation of the function to perform the transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date.
 *
 * @param xp Pole coordinate in x direction
 * @param yp Pole coordinate in y direction
 * @param[out] PoleMat Pole matrix
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\PoleMatrix.h"

Matrix PoleMatrix(double xp, double yp){
    Matrix rx = R_x(-yp);
    Matrix ry = R_y(-xp);
    
    return ry * rx;
}