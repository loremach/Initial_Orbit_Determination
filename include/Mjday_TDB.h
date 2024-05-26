// $Header$
//----------------------------------------------------------------------------------------
//                          Mjday_TDB
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/23
//
/*
 * @file Mjday_TDB.h
 * @brief Header file for computing the Modified Julian Date for barycentric dynamical time (TDB)
 *
 * @details This header file contains declarations related to computing the Modified Julian Date for barycentric dynamical time (TDB) from Modified Julian Date (TT).
 *
 * @param Mjd_TT Modified Julian Date (TT)
 * @return Mjd_TDB Modified Julian Date (TDB)
 *
 * @author Lorena Remacha Bordallo
*/
#ifndef _MJDAYTDB_
#define _MJDAYTDB_

#include <math.h>

/**
 * @brief Computes the Modified Julian Date for barycentric dynamical time (TDB).
 *
 * @param Mjd_TT Modified Julian Date (TT)
 * @return Mjd_TDB Modified Julian Date (TDB)
 */
double Mjday_TDB(double Mjd_TT);

#endif