// $Source$
//----------------------------------------------------------------------------------------
//                          hgibbs
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/21
//
/*
 * @file hgibbs.cpp
 * @brief Implements the Herrick-Gibbs approximation for orbit determination.
 *
 * @details This file contains the implementation of the function to implement the Herrick-Gibbs approximation
 * for orbit determination, and finds the middle velocity vector for the 3 given position vectors.
 *
 * @param r1 IJK position vector #1 [m]
 * @param r2 IJK position vector #2 [m]
 * @param r3 IJK position vector #3 [m]
 * @param Mjd1 Julian date of 1st sighting [days from 4713 BC]
 * @param Mjd2 Julian date of 2nd sighting [days from 4713 BC]
 * @param Mjd3 Julian date of 3rd sighting [days from 4713 BC]
 * @param[out] v2 IJK velocity vector for r2 [m/s]
 * @param[out] theta Angle between vectors [rad]
 * @param[out] error Flag indicating success ('ok',...)
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\hgibbs.h"

void hgibbs (Matrix r1, Matrix r2, Matrix r3, double Mjd1, double Mjd2, double Mjd3, Matrix & v2, double & theta, double &theta1, double & cop, string & error){
    error =  "ok";
    theta = 0.0;
    theta1= 0.0;
    double magr1 = norm( r1 );
    double magr2 = norm( r2 );
    double magr3 = norm( r3 );

    for (int i = 1; i <= 3; i++) {
        v2(i) = 0.0;
    }

    double tolangle= 0.01745329251994;
    double dt21= (Mjd2-Mjd1)*86400.0;
    double dt31= (Mjd3-Mjd1)*86400.0;
    double dt32= (Mjd3-Mjd2)*86400.0;

    Matrix p = cross( r2,r3 );
    Matrix pn = unit( p );
    Matrix r1n = unit( r1 );
    cop =  asin( dot( pn,r1n ) );

    if (fabs( dot(r1n,pn)) > 0.017452406 ) {
        error = "not coplanar";
    }

    theta  = angl( r1,r2 );
    theta1 = angl( r2,r3 );

    if ((theta > tolangle) | (theta1 > tolangle) )
        error= "angl > 1ø";


    double term1= -dt32*( 1.0/(dt21*dt31) + Const::GM_Earth/(12.0*magr1*magr1*magr1));
    double term2= (dt32-dt21)*( 1.0/(dt21*dt32) + Const::GM_Earth/(12.0*magr2*magr2*magr2));
    double term3=  dt21*( 1.0/(dt32*dt31) + Const::GM_Earth/(12.0*magr3*magr3*magr3));

    v2 =  r1*term1 + r2*term2 + r3*term3;
}