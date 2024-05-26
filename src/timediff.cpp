// $Source$
//----------------------------------------------------------------------------------------
//                          timediff
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/23
//
/*
 * @file timediff.cpp
 * @brief Time differences in seconds
 *
 * @details This file contains the implementation of the function to compute time differences in seconds.
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

#include "..\include\timediff.h"

void timediff(double UT1_UTC, double TAI_UTC, double& UT1_TAI, double& UTC_GPS, double& UT1_GPS, double& TT_UTC, double& GPS_UTC){
    double TT_TAI  = +32.184;          // TT-TAI time difference [s]

    double GPS_TAI = -19.0;            // GPS-TAI time difference [s]

    double TT_GPS  =  TT_TAI-GPS_TAI;  // TT-GPS time difference [s]

    double TAI_GPS = -GPS_TAI;         // TAI-GPS time difference [s]

    UT1_TAI = UT1_UTC-TAI_UTC;  // UT1-TAI time difference [s]

    double UTC_TAI = -TAI_UTC;         // UTC-TAI time difference [s]
    
    UTC_GPS = UTC_TAI-GPS_TAI;  // UTC_GPS time difference [s]

    UT1_GPS = UT1_TAI-GPS_TAI;  // UT1-GPS time difference [s]

    TT_UTC  = TT_TAI-UTC_TAI;   //  TT-UTC time difference [s]

    GPS_UTC = GPS_TAI-UTC_TAI;  // GPS-UTC time difference [s]
}