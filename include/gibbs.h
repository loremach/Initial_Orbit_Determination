// $Header$
//----------------------------------------------------------------------------------------
//                          gibbs
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/18
//
/*
 * @file gibbs.h
 * @brief Header file for the Gibbs method of orbit determination
 *
 * @details This header file contains the declaration of the function to perform the Gibbs method of orbit determination.
 *
 * @param r1 IJK position vector #1 [m]
 * @param r2 IJK position vector #2 [m]
 * @param r3 IJK position vector #3 [m]
 * @param[out] v2 IJK velocity vector for r2 [m/s]
 * @param[out] theta Angle between vectors [rad]
 * @param[out] error Flag indicating success ('ok',...)
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _GIBBS_
#define _GIBBS_

#include <math.h>
#include <string>
#include "..\include\Matrix.h"
#include "..\include\SAT_Const.h"
#include "..\include\unit.h"
#include "..\include\angl.h"

/**
 * @brief Performs the Gibbs method of orbit determination.
 *
 * @param r1 IJK position vector #1 [m]
 * @param r2 IJK position vector #2 [m]
 * @param r3 IJK position vector #3 [m]
 * @param[out] v2 IJK velocity vector for r2 [m/s]
 * @param[out] theta Angle between vectors [rad]
 * @param[out] error Flag indicating success ('ok',...)
 */
void gibbs(Matrix r1,Matrix r2,Matrix r3, Matrix & v2, double & theta, double & theta1, double & cop, string & error);

#endif