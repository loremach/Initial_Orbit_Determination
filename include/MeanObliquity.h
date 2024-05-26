// $Header$
//----------------------------------------------------------------------------------------
//                          MeanObliquity
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
 * @file MeanObliquity.h
 * @brief Header file for computing the mean obliquity of the ecliptic
 *
 * @details This header file contains declarations related to computing the mean obliquity of the ecliptic.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return MOblq Mean obliquity of the ecliptic [rad]
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _MEANOBLIQUITY_
#define _MEANOBLIQUITY_

#include <math.h>
#include "SAT_Const.h"

/**
 * @brief Computes the mean obliquity of the ecliptic.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return MOblq Mean obliquity of the ecliptic [rad]
 */
double MeanObliquity (double Mjd_TT);

#endif