// $Header$
//----------------------------------------------------------------------------------------
//                          VarEqn
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/21
//
/*
 * @file VarEqn.h
 * @brief Header file for computing the variational equations
 *
 * @details This header file contains declarations for computing the variational equations,
 * which involve the derivative of the state vector and the state transition matrix.
 *
 * @param x Time since epoch in [s]
 * @param yPhi (6+36)-dimensional vector comprising the state vector (y) and the state transition matrix (Phi) in column-wise storage order
 * @param[out] yPhip Derivative of yPhi
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef VAREQN
#define VAREQN

#include <math.h>
#include "..\include\Matrix.h"
#include "..\include\global.h"
#include "..\include\IERS.h"
#include "..\include\timediff.h"
#include "..\include\PrecMatrix.h"
#include "..\include\NutMatrix.h"
#include "..\include\PoleMatrix.h"
#include "..\include\GHAMatrix.h"
#include "..\include\AccelHarmonic.h"
#include "..\include\G_AccelHarmonic.h"

/**
 * @brief Computes the variational equations.
 *
 * @param x Time since epoch in [s]
 * @param yPhi (6+36)-dimensional vector comprising the state vector (y) and the state transition matrix (Phi) in column-wise storage order
 * @param[out] yPhip Derivative of yPhi
 */
Matrix VarEqn(double x, Matrix & yPhi);

#endif
