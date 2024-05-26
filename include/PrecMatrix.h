// $Header$
//----------------------------------------------------------------------------------------
//                          PrecMatrix
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/03
//
/*
 * @file PrecMatrix.h
 * @brief Header file for precession transformation of equatorial coordinates
 *
 * @details This header file contains declarations related to precession transformation of equatorial coordinates.
 *
 * @param Mjd_1 Epoch given (Modified Julian Date TT)
 * @param MjD_2 Epoch to precess to (Modified Julian Date TT)
 * @param[out] PrecMat Precession transformation matrix
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _PRECMATRIX_
#define _PRECMATRIX_

#include <math.h>
#include "SAT_const.h"
#include "..\include\R_z.h"
#include "..\include\R_y.h"
#include "..\include\Matrix.h"

/**
 * @brief Precession transformation of equatorial coordinates.
 *
 * @param Mjd_1 Epoch given (Modified Julian Date TT)
 * @param MjD_2 Epoch to precess to (Modified Julian Date TT)
 * @param[out] PrecMat Precession transformation matrix
 */
Matrix PrecMatrix(double Mjd_1, double Mjd_2);

#endif