// $Header$
//----------------------------------------------------------------------------------------
//                          timediff
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/23
//
/*
 * @file timediff.h
 * @brief Header file for computing time differences in seconds
 *
 * @details This header file contains declarations related to computing time differences in seconds.
 *
 * @param UT1_UTC Difference between UT1 and UTC [s]
 * @param TAI_UTC Difference between TAI and UTC [s]
 * @param[out] UT1_TAI Difference between UT1 and TAI [s]
 * @param[out] UTC_GPS Difference between UTC and GPS [s]
 * @param[out] UT1_GPS Difference between UT1 and GPS [s]
 * @param[out] TT_UTC Difference between Terrestrial Time (TT) and UTC [s]
 * @param[out] GPS_UTC Difference between GPS time and UTC [s]
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _TIMEDIFF_
#define _TIMEDIFF_

#include <math.h>

/**
 * @brief Computes various time differences.
 *
 * @param UT1_UTC Difference between UT1 and UTC [s]
 * @param TAI_UTC Difference between TAI and UTC [s]
 * @param[out] UT1_TAI Difference between UT1 and TAI [s]
 * @param[out] UTC_GPS Difference between UTC and GPS [s]
 * @param[out] UT1_GPS Difference between UT1 and GPS [s]
 * @param[out] TT_UTC Difference between Terrestrial Time (TT) and UTC [s]
 * @param[out] GPS_UTC Difference between GPS time and UTC [s]
 */
void timediff(double UT1_UTC, double TAI_UTC, double& UT1_TAI, double& UTC_GPS, double& UT1_GPS, double& TT_UTC, double& GPS_UTC);

#endif