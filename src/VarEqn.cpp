// $Source$
//----------------------------------------------------------------------------------------
//                          VarEqn
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/21
//
/*
 * @file VarEqn.cpp
 * @brief Header file for computing the variational equations
 *
 * @details This file contains the implementation for computing the variational equations,
 * which involve the derivative of the state vector and the state transition matrix.
 *
 * @param x Time since epoch in [s]
 * @param yPhi (6+36)-dimensional vector comprising the state vector (y) and the state transition matrix (Phi) in column-wise storage order
 * @param[out] yPhip Derivative of yPhi
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\VarEqn.h"

Matrix VarEqn(double x, Matrix & yPhi){
    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;
    IERS(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, *Global::eopdata, Global::auxparam.Mjd_UTC, 'l');

    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    double Mjd_UT1 = Global::auxparam.Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    // Transformation matrix
    Matrix P = PrecMatrix(Const::MJD_J2000,Global::auxparam.Mjd_TT + x/86400.0);
    Matrix N = NutMatrix(Global::auxparam.Mjd_TT + x/86400);
    Matrix T = N * P;

    Matrix pole = PoleMatrix(x_pole,y_pole);
    Matrix gha = GHAMatrix(Mjd_UT1);
    Matrix E = pole * gha * T;

    // State vector components
    Matrix r(3), v(3);
    for (int i = 1; i <= 3; i++){
        r(i) = yPhi(i);
        v(i) = yPhi(i+3);
    }

    Matrix Phi = zeros(6,6);

    // State transition matrix
    Matrix aux(6);
    for (int j = 1; j <= 6; j++){
        for(int i =1; i <= 6; i++){
            aux(i) = yPhi(6*j+i);
        }
        Phi = assign_column(Phi, j, transpose(aux));
    }

    // Acceleration and gradient
    Matrix a = AccelHarmonic ( r, E, Global::auxparam.n, Global::auxparam.m );
    Matrix G = G_AccelHarmonic ( r, E, Global::auxparam.n, Global::auxparam.m );
    // Time derivative of state transition matrix
    Matrix yPhip = zeros(1,42);
    Matrix df_dy = zeros(6, 6);

    for(int i = 1; i <= 3; i++){
        for(int j = 1; j <= 3; j++){
            df_dy(i,j) = 0.0;
            df_dy(i+3,j) = G(i,j);
            if(i==j){
                df_dy(i,j+3) = 1.0;
            }else{
                df_dy(i,j+3) = 0.0;
            }
            df_dy(i+3,j+3) = 0.0;
        }
    }

    Matrix Phip = zeros(6,6);
    Phip = df_dy*Phi;

    // Derivative of combined state vector and state transition matrix
    for(int i = 1; i <= 3; i++){
        yPhip(i)   = v(i);                 // dr/dt(i)
        yPhip(i+3) = a(i);
    }

    for(int i = 1; i <= 6; i++){
        for(int j = 1; j <= 6; j++) {
            yPhip(6*j+i) = Phip(i,j);     // dPhi/dt(i,j)
        }
    }

    return yPhip;
}

