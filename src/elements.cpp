// $Source$
//----------------------------------------------------------------------------------------
//                          elements
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/05
//
/*
 * @file elements.cpp
 * @brief Computes the osculating Keplerian elements from the satellite state vector for elliptic orbits
 *
 * @details This file contains the declaration of the function to compute the osculating Keplerian elements
 * from the satellite state vector for elliptic orbits. The function cannot be used with state vectors
 * describing a circular or non-inclined orbit.
 *
 * @param y State vector (x, y, z, vx, vy, vz)
 * @param p Semilatus rectum in meters
 * @param a Semimajor axis
 * @param e Eccentricity
 * @param i Inclination in radians
 * @param Omega Longitude of the ascending node in radians
 * @param omega Argument of pericenter in radians
 * @param M Mean anomaly in radians
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\elements.h"
#include "..\include\utils.h"

void elements(double & p, double & a, double & e, double & i, double & Omega, double & omega, double & M, Matrix y){ 
    Matrix r(3);    r(1) = y(1);    r(2) = y(2);    r(3) = y(3);    // Position
    Matrix v(3);    v(1) = y(4);    v(2) = y(5);    v(3) = y(6);    // Velocity                                     

    Matrix h = cross(r,v);                                    // Areal velocity
    double magh = norm(h);
    p = magh*magh/Const::GM_Earth;
    double H = norm(h);
    
    Omega = atan2 ( h(1), -h(2) );  
    Omega = custom_mod(Omega, Const::pi2);          
    i     = atan2 ( sqrt(h(1)*h(1)+h(2)*h(2)), h(3) ); // Inclination        
    double u     = atan2 ( r(3)*H, -r(1)*h(2)+r(2)*h(1) );    // Arg. of latitude   
    
    double R  = norm(r);                                      // Distance           
    
    a = 1/(2/R-dot(v,v)/Const::GM_Earth);               // Semi-major axis    
    
    double eCosE = 1-R/a;                                     // e*cos(E)           
    double eSinE = dot(r,v)/sqrt(Const::GM_Earth*a);           // e*sin(E)           
    
    double e2 = eCosE*eCosE +eSinE*eSinE;
    e  = sqrt(e2);                                     // Eccentricity 
    double E  = atan2(eSinE,eCosE);                           // Eccentric anomaly  
    
    M  = custom_mod(E-eSinE, Const::pi2);                             // Mean anomaly
    
    double nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);          // True anomaly
    
    omega = custom_mod(u-nu,Const::pi2);                             // Arg. of perihelion 
}

