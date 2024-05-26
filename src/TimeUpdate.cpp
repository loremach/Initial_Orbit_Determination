// $Source$
//----------------------------------------------------------------------------------------
//                          TimeUpdate
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/03
//
/*
 * @file TimeUpdate.cpp
 *
 * @brief File for time update function
 *
 * @param P Covariance matrix
 * @param Phi State transition matrix
 * @param Qdt Process noise covariance matrix multiplied by time step
 * @return Updated covariance matrix
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\TimeUpdate.h"

Matrix TimeUpdate(Matrix P, Matrix Phi, double Qdt){
    return Phi*P*transpose(Phi) + Qdt;
}