// $Header$
//----------------------------------------------------------------------------------------
//                          MeasUpdate
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/04
//
/*
 * @file MeasUpdate.h
 * @brief Header file for performing measurement update in a Kalman filter.
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

#ifndef _MEASUPDATE_
#define _MEASUPDATE_

#include <math.h>
#include "matrix.h"

/**
 * @brief Function to perform measurement update in a Kalman filter.
 * @param K Kalman gain matrix (output).
 * @param x State vector (input/output).
 * @param P State covariance matrix (input/output).
 * @param z Measurement value.
 * @param g Predicted measurement value.
 * @param s Measurement standard deviation.
 * @param G Measurement sensitivity matrix.
 * @param n Size of the state vector.
 */
void MeasUpdate(Matrix & K, Matrix & x, Matrix & P, double z, double g, double s, Matrix G, double n);

#endif