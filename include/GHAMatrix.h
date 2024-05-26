// $Header$
//----------------------------------------------------------------------------------------
//                          GHAMatrix
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/05
//
/*
 * @file GHAMatrix.h
 * @brief Header file for Transformation from true equator and equinox to Earth equator and Greenwich meridian system
 *
 * @details This header file contains the declaration of the function to compute the transformation from true equator and equinox
 * to Earth equator and Greenwich meridian system.
 *
 * @param Mjd_UT1 Modified Julian Date UT1
 * @param[out] GHAmat Greenwich Hour Angle matrix
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _GHAMATRIX_
#define _GHAMATRIX_

#include <math.h>
#include "matrix.h"
#include "R_z.h"
#include "gast.h"

/**
 * @brief Transformation from true equator and equinox to Earth equator and Greenwich meridian system.
 *
 * @param Mjd_UT1 Modified Julian Date UT1
 * @param[out] GHAmat Greenwich Hour Angle matrix
 */
Matrix GHAMatrix (double Mjd_UT1);

#endif