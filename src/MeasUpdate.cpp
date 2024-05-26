// $Source$
//----------------------------------------------------------------------------------------
//                          MeasUpdate
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/04
//
/*
 * @file MeasUpdate.cpp
 * @brief Perform measurement update in a Kalman filter.
 *
 * @param K Kalman gain matrix (output).
 * @param x State vector (input/output).
 * @param P State covariance matrix (input/output).
 * @param z Measurement value.
 * @param g Predicted measurement value.
 * @param s Measurement standard deviation.
 * @param G Measurement sensitivity matrix.
 * @param n Size of the state vector.
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\MeasUpdate.h"

void MeasUpdate(Matrix & K, Matrix & x, Matrix & P, double z, double g, double s, Matrix G, double n){
    double m=1;

    Matrix Inv_W = zeros(m,m);
    Inv_W(1, 1) = s*s;    // Inverse weight (measurement covariance)

    // Kalman gain
    K = transpose(P*transpose(G)*inv(Inv_W+G*P*transpose(G)));

    // State update
    x = x + K*(z-g);

    // Covariance update
    P = (eye(n)-transpose(K)*G)*P;
}