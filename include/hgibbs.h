// $Header$
//----------------------------------------------------------------------------------------
//                          hgibbs
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/21
//
/*
 * @file hgibbs.h
 * @brief Header file for implementing the Herrick-Gibbs approximation for orbit determination
 *
 * @details This header file contains the declaration of the function to implement the Herrick-Gibbs approximation
 * for orbit determination, and finds the middle velocity vector for the 3 given position vectors.
 *
 * @param r1 IJK position vector #1 [m]
 * @param r2 IJK position vector #2 [m]
 * @param r3 IJK position vector #3 [m]
 * @param Mjd1 Julian date of 1st sighting [days from 4713 BC]
 * @param Mjd2 Julian date of 2nd sighting [days from 4713 BC]
 * @param Mjd3 Julian date of 3rd sighting [days from 4713 BC]
 * @param[out] v2 IJK velocity vector for r2 [m/s]
 * @param[out] theta Angle between vectors [rad]
 * @param[out] error Flag indicating success ('ok',...)
 *
 * @author Lorena Remacha Bordallo
*/
#ifndef HGIBBS
#define HGIBBS

#include <math.h>
#include <string>
#include "..\include\Matrix.h"
#include "..\include\SAT_Const.h"
#include "..\include\unit.h"
#include "..\include\angl.h"

/**
 * @brief Implements the Herrick-Gibbs approximation for orbit determination.
 *
 * @param r1 IJK position vector #1 [m]
 * @param r2 IJK position vector #2 [m]
 * @param r3 IJK position vector #3 [m]
 * @param Mjd1 Julian date of 1st sighting [days from 4713 BC]
 * @param Mjd2 Julian date of 2nd sighting [days from 4713 BC]
 * @param Mjd3 Julian date of 3rd sighting [days from 4713 BC]
 * @param[out] v2 IJK velocity vector for r2 [m/s]
 * @param[out] theta Angle between vectors [rad]
 * @param[out] error Flag indicating success ('ok',...)
 */
void hgibbs (Matrix r1, Matrix r2, Matrix r3, double Mjd1, double Mjd2, double Mjd3, Matrix & v2, double & theta, double &theta1, double & cop, string & error);
#endif
