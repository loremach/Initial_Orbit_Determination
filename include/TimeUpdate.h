// $Header$
//----------------------------------------------------------------------------------------
//                          TimeUpdate
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/03
//
/*
 * @file TimeUpdate.h
 *
 * @brief Header file for time update function
 *
 * @details This header file contains declarations related to the time update function.
 *
 * @param P Covariance matrix
 * @param Phi State transition matrix
 * @param Qdt Process noise covariance matrix multiplied by time step
 * @return Updated covariance matrix
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include <math.h>
#include "..\include\Matrix.h"

/**
 * @brief Performs time update of the covariance matrix.
 *
 * @param P Covariance matrix
 * @param Phi State transition matrix
 * @param Qdt Process noise covariance matrix multiplied by time step
 * @return Updated covariance matrix
 */
Matrix TimeUpdate(Matrix P, Matrix Phi, double Qdt=0.0);

#endif