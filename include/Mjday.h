// $Header$
//----------------------------------------------------------------------------------------
//                          Mjday
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/23
//
/*
 * @file Mjday.h
 * @brief Header file for computing the Modified Julian Date (Mjd)
 *
 * @details This header file contains declarations related to computing the Modified Julian Date (Mjd) from year, month, day, hour, minute, and second.
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


#ifndef _MJDAY_
#define _MJDAY_

#include <math.h>

/**
 * @brief Computes the Modified Julian Date (Mjd) from year, month, day, hour, minute, and second.
 *
 * @param year Year
 * @param mon Month
 * @param day Day
 * @param hr Universal time hour
 * @param min Universal time minute
 * @param sec Universal time second
 * @return Mjd Modified Julian Date
 */
double Mjday(int yr, int mon, int day, int hr=0, int min=0, int sec=0);

#endif