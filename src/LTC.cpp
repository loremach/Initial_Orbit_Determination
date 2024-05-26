// $Source$
//----------------------------------------------------------------------------------------
//                          LTC
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/04
//
/*
 * @file LTC.cpp
* @brief Transformation from Greenwich meridian system to local tangent coordinates
 *
 * @details This file contains the implementation of the function to perform the transformation from Greenwich meridian system to local tangent coordinates.
 *
 * @param lon Geodetic East longitude [rad]
 * @param lat Geodetic latitude [rad]
 * @param[out] M Rotation matrix from the Earth equator and Greenwich meridian to the local tangent (East-North-Zenith) coordinate system
 *
 * @author Lorena Remacha Bordallo
*/
#include "..\include\LTC.h"

Matrix LTC(double lon, double lat){
    Matrix r1 = R_y(-1.0*lat);
    Matrix r2 = R_z(lon);
    Matrix M = r1*r2;
    double Aux;
    for (int j=1; j<=3; j++){
        Aux=M(1,j); M(1,j)=M(2,j); M(2,j)=M(3,j); M(3,j)= Aux;
    }
    return M;
}