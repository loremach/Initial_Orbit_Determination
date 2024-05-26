// $Source$
//----------------------------------------------------------------------------------------
//                          R_x
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/11
//
/*
 * @file R_x.cpp
 * @brief Computes the vector result of rotation
 *
 * @details This file contains the implementation of the function to compute the vector result of rotation.
 *
 * @param angle Angle of rotation [rad]
 * @param[out] rotmat Vector result
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\R_x.h"

Matrix R_x(double angle){
    double C, S;

    C = cos(angle);
    S = sin(angle);
    Matrix rotmat = zeros(3,3);

    rotmat(1,1) = 1.0;  rotmat(1,2) =    0.0;  rotmat(1,3) = 0.0;
    rotmat(2,1) = 0.0;  rotmat(2,2) =      C;  rotmat(2,3) =   S;
    rotmat(3,1) = 0.0;  rotmat(3,2) = -1.0*S;  rotmat(3,3) =   C;

    return rotmat;
}