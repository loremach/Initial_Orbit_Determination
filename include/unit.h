// $Header$
//----------------------------------------------------------------------------------------
//                          unit
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/23
//
/*
 * @file unit.h
 * @brief Header file for computing unit vectors
 *
 * @details This header file contains declarations related to computing unit vectors.
 *
 * @param vec Input vector
 * @param[out] outvec Unit vector
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _UNIT_
#define _UNIT_

#include "matrix.h"
#include <math.h>

/**
 * @brief Computes the unit vector given the original vector.
 *
 * @param vec Input vector
 * @param[out] outvec Unit vector
 */
Matrix unit(Matrix vec);

#endif