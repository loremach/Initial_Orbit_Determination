// $Source$
//----------------------------------------------------------------------------------------
//                          GHAMatrix
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/05
//
/*
 * @file GHAMatrix.cpp
 * @brief Transformation from true equator and equinox to Earth equator and Greenwich meridian system
 *
 * @details This file contains the implementation of the function to compute the transformation from true equator and equinox
 * to Earth equator and Greenwich meridian system.
 *
 * @param Mjd_UT1 Modified Julian Date UT1
 * @param[out] GHAmat Greenwich Hour Angle matrix
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\GHAMatrix.h"

Matrix GHAMatrix (double Mjd_UT1){
    return R_z( gast(Mjd_UT1) );
}