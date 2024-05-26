// $Header$
//----------------------------------------------------------------------------------------
//                          global
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/18
//
/*
 * @file global.h
 * @brief Header file for global data initialization
 *
 * @details This file contains the implementation of functions to initialize global data
 * required by the application, such as Earth Orientation Parameters (EOP), gravity field coefficients, etc.
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _GLOBAL_
#define _GLOBAL_

#include "matrix.h"
#include "Mjday.h"
#include "SAT_Const.h"
#include <cstdio>
#include <cstdlib>
#include <string.h>

typedef struct{
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
} Param;

class Global{
public:
    static Param auxparam;
    static Matrix *eopdata;
    static Matrix *Cnm;
    static Matrix *Snm;
    static Matrix *PC;
    static Matrix *obs;

    /**
     * @brief Initializes the Earth Orientation Parameters (EOP) data.
     *
     * @param c Number of columns in the EOP data matrix.
     */
    static void eop19620101(int c);

    /**
     * @brief Initializes the gravity field coefficients from GGM03S model.
     */
    static void GGM03S();

    /**
     * @brief Initializes the coefficients for DE430 planetary ephemeris model.
     *
     * @param f Number of rows in the coefficient matrix.
     * @param c Number of columns in the coefficient matrix.
     */
    static void DE430Coeff(int f, int c);

    /**
     * @brief Initializes the GEOS-3 satellite observations.
     *
     * @param f Number of observations.
     */
    static void GEOS3(int f);
};

#endif