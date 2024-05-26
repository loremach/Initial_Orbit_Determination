// $Header$
//----------------------------------------------------------------------------------------
//                          PoleMatrix
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/03
//
/*
 * @file PoleMatrix.h
 * @brief Header file for transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date
 *
 * @details This header file contains declarations related to transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date.
 *
 * @param xp Pole coordinate in x direction
 * @param yp Pole coordinate in y direction
 * @param[out] PoleMat Pole matrix
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _POLEMATRIX_
#define _POLEMATRIX_

#include <math.h>
#include "..\include\R_x.h"
#include "..\include\R_y.h"
#include "..\include\Matrix.h"

/**
 * @brief Transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date.
 *
 * @param xp Pole coordinate in x direction
 * @param yp Pole coordinate in y direction
 * @param[out] PoleMat Pole matrix
 */
Matrix PoleMatrix(double xp, double yp);

#endif