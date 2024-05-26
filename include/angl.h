// $Header$
//----------------------------------------------------------------------------------------
//                          angl
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/04
//
/*
* @file angl.h
* @brief Header file for computing the angle between two vectors
 *
 * @details This header file contains the implementation of the function to compute the angle between two vectors.
 *
 * @param vec1 Vector 1
 * @param vec2 Vector 2
 * @return theta Angle between the two vectors in the range -pi to pi
 *
 * @author Lorena Remacha Bordallo
*/
#ifndef _ANGL_
#define _ANGL_

#include <math.h>
#include "Matrix.h"

/**
 * @brief Computes the angle between two vectors.
 *
 * @param vec1 Vector 1
 * @param vec2 Vector 2
 * @return theta Angle between the two vectors in the range -pi to pi
 */
double angl (Matrix vec1, Matrix vec2);

#endif