// $Header$
//----------------------------------------------------------------------------------------
//                          NutMatrix
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/02
//
/*
 * @file NutMatrix.h
 * @brief Header file for transformation from mean to true equator and equinox
 *
 * @details This header file contains declarations related to transformation from mean to true equator and equinox.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @param[out] NutMat Nutation matrix
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _NUTMATRIX_
#define _NUTMATRIX_

#include <math.h>
#include "..\include\NutAngles.h"
#include "..\include\MeanObliquity.h"
#include "..\include\R_x.h"
#include "..\include\R_y.h"
#include "..\include\R_z.h"
#include "..\include\Matrix.h"

/**
 * @brief Transformation from mean to true equator and equinox.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @param[out] NutMat Nutation matrix
 */
Matrix NutMatrix (double Mjd_TT);

#endif