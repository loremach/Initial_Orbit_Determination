// $Source$
//----------------------------------------------------------------------------------------
//                          Cheb3D
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
 * @file Cheb3D.cpp
 * @brief Computes Chebyshev approximation of 3-dimensional vectors
 *
 * @details This file contains the declaration of the function for Chebyshev approximation
 * of 3-dimensional vectors.
 *
 * @param N Number of coefficients
 * @param Ta Begin interval
 * @param Tb End interval
 * @param Cx Coefficients of Chebyshev polynomial (x-coordinate)
 * @param Cy Coefficients of Chebyshev polynomial (y-coordinate)
 * @param Cz Coefficients of Chebyshev polynomial (z-coordinate)
 * @param t The evaluation point within the interval [Ta, Tb]
 * @param result The resulting 3-dimensional vector approximation
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\Cheb3D.h"

Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix Cx, Matrix Cy, Matrix Cz){
    // Check validity
    if ( (t<Ta) || (Tb<t) ){
        cout<<"ERROR: Time out of range in Cheb3D::Value\n";
        exit(EXIT_FAILURE);
    }

    // Clenshaw algorithm
    double tau = (2*t-Ta-Tb)/(Tb-Ta);  

    Matrix f1 = zeros(1, 3);
    Matrix f2 = zeros(1, 3);

    Matrix old_f1(3);
    Matrix C(3); 

    for (int i = N; i >= 2; --i) {
        old_f1 = f1;
        C(1) = Cx(i);     C(2) = Cy(i);    C(3) = Cz(i);
        f1 = (f1*tau)*2 - f2 + C;
        f2 = old_f1;
    }
    C(1) = Cx(1);     C(2) = Cy(1);    C(3) = Cz(1);
    return f1*tau-f2+C;
}