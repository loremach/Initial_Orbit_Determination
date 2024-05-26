// $Source$
//----------------------------------------------------------------------------------------
//                          G_AccelHarmonic
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/08
//
/*
 * @file G_AccelHarmonic.cpp
 * @brief Computes the gradient of the Earth's harmonic gravity field
 *
 * @details This file contains the implementation of the function to compute the gradient of the Earth's harmonic gravity field.
 *
 * @param r Satellite position vector in the true-of-date system
 * @param U Transformation matrix to body-fixed system
 * @param n Gravity model degree
 * @param m Gravity model order
 * @return G Gradient (G = da/dr) in the true-of-date system
 *
 * @author Lorena Remacha Bordallo
*/


#include "..\include\G_AccelHarmonic.h"

Matrix G_AccelHarmonic(Matrix r, Matrix U, int n_max, int m_max){
    double d = 1.0;   // Position increment [m]

    Matrix G = zeros(3, 3);
    Matrix dr = zeros(1, 3);

    Matrix da(1, 3);
    Matrix aux1(1, 3);
    Matrix aux2(1, 3);

    // Gradient
    for (int i=1; i<=3; i++){
        // Set offset in i-th component of the position vector
        dr = zeros(1, 3);
        dr(i) = d;
        
        aux1 = AccelHarmonic ( r+dr/2,U, n_max, m_max );
        aux2 = AccelHarmonic ( r-dr/2,U, n_max, m_max );

        da = aux1 - aux2;
        
        G = assign_column(G, i, transpose(da/d));

    }

    return G;
}