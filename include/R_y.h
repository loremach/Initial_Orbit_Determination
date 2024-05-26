// $Header$
//----------------------------------------------------------------------------------------
//                          R_y
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/14
//
/*
 * @file R_y.h
 * @brief Header file for computing the vector result of rotation
 *
 * @details This header file contains declarations related to computing the vector result of rotation.
 *
 * @param angle Angle of rotation [rad]
 * @param[out] rotmat Vector result
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _RY_
#define _RY_

#include "matrix.h"
#include <math.h>

/**
 * @brief Computes the vector result of rotation.
 *
 * @param angle Angle of rotation [rad]
 * @param[out] rotmat Vector result
 */
Matrix R_y(double angle);

#endif