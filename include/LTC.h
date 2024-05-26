// $Header$
//----------------------------------------------------------------------------------------
//                          LTC
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/04
//
/*
 * @file LTC.h
 * @brief Header file for transformation from Greenwich meridian system to local tangent coordinates
 *
 * @details This header file contains declarations related to transformation from Greenwich meridian system to local tangent coordinates.
 *
 * @param lon Geodetic East longitude [rad]
 * @param lat Geodetic latitude [rad]
 * @param[out] M Rotation matrix from the Earth equator and Greenwich meridian to the local tangent (East-North-Zenith) coordinate system
 *
 * @author Lorena Remacha Bordallo
*/
#ifndef _LTC_
#define _LTC_

#include <math.h>
#include "..\include\Matrix.h"
#include "..\include\R_z.h"
#include "..\include\R_y.h"

/**
 * @brief Transformation from Greenwich meridian system to local tangent coordinates.
 *
 * @param lon Geodetic East longitude [rad]
 * @param lat Geodetic latitude [rad]
 * @param[out] M Rotation matrix from the Earth equator and Greenwich meridian to the local tangent (East-North-Zenith) coordinate system
 */
Matrix LTC(double lon, double lat);

#endif