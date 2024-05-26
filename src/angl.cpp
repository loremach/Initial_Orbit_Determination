// $Source$
//----------------------------------------------------------------------------------------
//                          angl
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/04
//
/*
 * @file angl.cpp
 * @brief Computes the angle between two vectors
 *
 * @details This file contains the implementation of the function to compute the angle between two vectors.
 *
 * @param vec1 Vector 1
 * @param vec2 Vector 2
 * @return theta Angle between the two vectors in the range -pi to pi
 *
 * @author Lorena Remacha Bordallo
*/
#include "..\include\angl.h"

double angl(Matrix vec1, Matrix vec2)
{
    double small = 0.00000001;
    double undefined = 999999.1;

    double magv1 = norm(vec1);
    double magv2 = norm(vec2);
    double temp;
    double theta;

    if ((magv1 * magv2) > pow(small,2))
    {
        temp = dot(vec1, vec2) / (magv1 * magv2);
        if (abs(temp) > 1.0)
        {
            if (temp >= 0)
                temp = temp * 1.0;
            else
                temp = -temp * 1.0;
        }
        theta = acos(temp);
    }
    else
    {
        theta = undefined;
    }
    return theta;
}