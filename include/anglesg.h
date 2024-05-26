// $Header$
//----------------------------------------------------------------------------------------
//                          anglesg
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/04
//
/*
 * @file anglesg.h
 * @brief Header file for solving the problem of orbit determination using three optical sightings
 *
 * @param az1 Azimuth at t1 in radians
 * @param az2 Azimuth at t2 in radians
 * @param az3 Azimuth at t3 in radians
 * @param el1 Elevation at t1 in radians
 * @param el2 Elevation at t2 in radians
 * @param el3 Elevation at t3 in radians
 * @param Mjd1 Modified Julian Date of t1
 * @param Mjd2 Modified Julian Date of t2
 * @param Mjd3 Modified Julian Date of t3
 * @param Rs1 IJK site1 position vector in meters
 * @param Rs2 IJK site2 position vector in meters
 * @param Rs3 IJK site3 position vector in meters
 * @param r IJK position vector at t2 in meters
 * @param v IJK velocity vector at t2 in meters per second
 *
 * @author Lorena Remacha Bordallo
*/
#ifndef ANGLESG
#define ANGLESG

#include "Matrix.h"
#include "Geodetic.h"
#include "LTC.h"
#include "IERS.h"
#include "Global.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "SAT_const.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GhaMatrix.h"
#include "rpoly.h"
#include "gibbs.h"
#include "hgibbs.h"
#include "elements.h"

/**
 * @brief Solves the problem of orbit determination using three optical sightings.
 *
 * @param az1 Azimuth at t1 in radians
 * @param az2 Azimuth at t2 in radians
 * @param az3 Azimuth at t3 in radians
 * @param el1 Elevation at t1 in radians
 * @param el2 Elevation at t2 in radians
 * @param el3 Elevation at t3 in radians
 * @param Mjd1 Modified Julian Date of t1
 * @param Mjd2 Modified Julian Date of t2
 * @param Mjd3 Modified Julian Date of t3
 * @param Rs1 IJK site1 position vector in meters
 * @param Rs2 IJK site2 position vector in meters
 * @param Rs3 IJK site3 position vector in meters
 * @param r IJK position vector at t2 in meters
 * @param v IJK velocity vector at t2 in meters per second
 */
void anglesg (double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix Rs1, Matrix Rs2, Matrix Rs3, Matrix & r2, Matrix & v2);

#endif
