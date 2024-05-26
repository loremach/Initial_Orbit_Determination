// $Source$
//----------------------------------------------------------------------------------------
//                          unit
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/23
//
/*
 * @file unit.cpp
 * @brief Computes unit vectors
 *
 * @details This file contains the implementation of the function to compute unit vectors.
 *
 * @param vec Input vector
 * @param[out] outvec Unit vector
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\unit.h"

Matrix unit(Matrix vec){
    double small = 0.000001;
    double magv = norm(vec);
    Matrix outvec(3);

    if (magv > small){
        for (int i=1; i<=3; i++){
            outvec(i)= vec(i)/magv;
        }
    }else{
        for (int i=1; i<=3; i++){
            outvec(i)= 0.0;
        }
    }

    return outvec;
}