// $Source$
//----------------------------------------------------------------------------------------
//                          Mjday
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/23
//
/*
 * @file Mjday.cpp
 * @brief Computes the Modified Julian Date (Mjd)
 *
 * @details This file contains the implementation of the function to compute the Modified Julian Date (Mjd) from year, month, day, hour, minute, and second.
 *
 * @param year Year
 * @param mon Month
 * @param day Day
 * @param hr Universal time hour
 * @param min Universal time minute
 * @param sec Universal time second
 * @return Mjd Modified Julian Date
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\Mjday.h"

double Mjday(int yr, int mon, int day, int hr, int min, int sec){

    double jd = 367.0 * (double)yr
            - floor( (7.0 * ((double)yr + floor( ((double)mon + 9.0) / 12.0) ) ) * 0.25 )
            + floor( 275.0 * (double)mon / 9.0 )
            + (double)day + 1721013.5
            + ( ((double)sec/60.0 + (double)min ) / 60.0 + (double)hr ) / 24.0;

    return jd-2400000.5;
}