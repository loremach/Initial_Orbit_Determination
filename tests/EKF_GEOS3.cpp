// $Source$
//----------------------------------------------------------------------------------------
//                          EKF_GEOS3
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/23
//
/* @file initial_orbit_determination.c
 *  @brief Initial Orbit Determination using Gauss and Extended Kalman Filter methods
 *
 * References:
 * - O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and Applications", Springer Verlag, Heidelberg, 2000.
 * - D. Vallado, "Fundamentals of Astrodynamics and Applications", 4th Edition, 2013.
 * - G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.
 * - Meysam Mahooti
 *
 *  @author Lorena Remacha Bordallo
 *  @bug No known bugs
*/

#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>

#include "..\include\global.h"
#include "..\include\Matrix.h"
#include "..\include\SAT_const.h"
#include "..\include\Position.h"
#include "..\include\anglesg.h"
#include "..\include\Accel.h"
#include "..\include\DEInteg.hpp"
#include "..\include\LTC.h"
#include "..\include\IERS.h"
#include "..\include\timediff.h"
#include "..\include\VarEqn.h"
#include "..\include\gmst.h"
#include "..\include\TimeUpdate.h"
#include "..\include\AzElPa.h"
#include "..\include\MeasUpdate.h"

using namespace std;

int main(){
    Global::eop19620101(21413);
    Global::DE430Coeff(2285, 1020);
    Global::GEOS3(46);
    Global::GGM03S();
    int nobs = 46;

    double sigma_range = 92.5;          // [m]
    double sigma_az = 0.0224*Const::Rad; // [rad]
    double sigma_el = 0.0139*Const::Rad; // [rad]

    // Kaena Point station
    double lat = Const::Rad*21.5748;     // [rad]
    double lon = Const::Rad*(-158.2706); // [rad]
    double alt = 300.20;                // [m]

    Matrix Rs = Position(lon, lat, alt);

    double Mjd1 = (*Global::obs)(1,1);
    double Mjd2 = (*Global::obs)(9,1);
    double Mjd3 = (*Global::obs)(18,1);

    Matrix r2(3);
    Matrix v2(3);

    anglesg((*Global::obs)(1,2),(*Global::obs)(9,2),
            (*Global::obs)(18,2),(*Global::obs)(1,3),
            (*Global::obs)(9,3),(*Global::obs)(18,3),
            Mjd1,Mjd2,Mjd3,Rs,Rs,Rs,r2,v2);

    Matrix Y0_apr(6);
    for(int i = 1; i <= 3; i++){
        Y0_apr(i) = r2(i);
        Y0_apr(i+3) = v2(i);
    }

    double Mjd0 = Mjday(1995,1,29,02,38,0);

    double Mjd_UTC = (*Global::obs)(9, 1);
    Global::auxparam.Mjd_UTC = Mjd_UTC;
    Global::auxparam.n = 20;
    Global::auxparam.m = 20;
    Global::auxparam.sun = 1;
    Global::auxparam.moon = 1;
    Global::auxparam.planets = 1;

    int n_eqn  = 6;

    Matrix Y = DEInteg(Accel, n_eqn, Y0_apr, 0, -((*Global::obs)(9,1)-Mjd0)*86400.0, 1e-13, 1e-6);

    Matrix P = zeros(6, 6);

    for(int i = 1; i <= 3; i++){
        P(i,i)=1e8;
    }
    for(int i = 4; i <= 6; i++){
        P(i,i)=1e3;
    }

    Matrix LT = LTC(lon,lat);

    Matrix yPhi = zeros(1, 42);
    Matrix Phi  = zeros(6, 6);

    // Measurement loop
    double t = 0;

    double t_old;
    Matrix Y_old(6);

    Matrix eop = *Global::eopdata;
    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;

    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;

    double Mjd_TT, Mjd_UT1;
    double theta;
    double Azim;  double Elev;
    double Dist;

    Matrix r(3);
    Matrix U(3);
    Matrix s(3);
    Matrix dAds(3);
    Matrix dEds(3);
    Matrix aux1(3);
    Matrix aux2(3);
    Matrix aux3(3);
    Matrix aux4(3);
    Matrix dAdY(6);
    Matrix dEdY(6);
    Matrix K(6);
    Matrix dDds(3);
    Matrix dDdY(6);

    for (int i = 1; i <= nobs; i++){
        // Previous step
        t_old = t;
        Y_old = Y;
        // Time increment and propagation
        Mjd_UTC = (*Global::obs)(i,1); // Modified Julian Date
        t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]

        IERS(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, eop, Mjd_UTC, 'l');

        timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

        Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
        Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
        Global::auxparam.Mjd_UTC = Mjd_UTC;
        Global::auxparam.Mjd_TT = Mjd_TT;

        for (int ii = 1; ii <= 6; ii++) {
            yPhi(ii) = Y_old(ii);
            for (int j = 1; j <= 6; j++) {
                if (ii == j) {
                    yPhi(6 * j + ii) = 1;
                } else {
                    yPhi(6 * j + ii) = 0;
                }
            }
        }
        n_eqn = 42;
        yPhi = DEInteg (VarEqn, n_eqn, yPhi, 0, t-t_old, 1e-13, 1e-6);

        // Extract state transition matrices
        for(int j = 1; j<=6; j++){
            for(int k = 1; k <=6; k++){
                Phi(k,j) = yPhi(6*j+k);
            }
        }
        n_eqn = 6;

        Y = DEInteg (Accel, n_eqn, Y_old, 0, t-t_old, 1e-13, 1e-6);

        // Topocentric coordinates
        theta = gmst(Mjd_UT1);                    // Earth rotation
        U = R_z(theta);
        for (int k = 1; k <=3; k++){
            r(k) = Y(k);
        }

        s = transpose(LT*(U*transpose(r)-transpose(Rs)));                          // Topocentric position [m]

        // Time update
        P = TimeUpdate(P, Phi);

        // Azimuth and partials

        AzElPa(s, Azim, Elev, dAds, dEds);// Azimuth, Elevation
        aux1 = dAds*LT*U;
        aux2 = zeros(1,3);

        for( int k = 1; k <= 3; k++){
            dAdY(k) = aux1(k);
            dAdY(k+3) = aux2(k);
        }

        // Measurement update
        MeasUpdate(K, Y, P, (*Global::obs)(i,2), Azim, sigma_az, dAdY, 6);

        // Range and partials
        for (int k = 1; k <=3; k++){
            r(k) = Y(k);
        }

        s = transpose(LT*(U*transpose(r)-transpose(Rs)));                          // Topocentric position [m]

        AzElPa(s, Azim, Elev, dAds, dEds); // Azimuth, Elevation
        aux1 = dEds*LT*U;
        aux2 = zeros(1,3);
        for( int k = 1; k <= 3; k++){
            dEdY(k) = aux1(k);
            dEdY(k+3) = aux2(k);
        }

        // Measurement update
        MeasUpdate(K, Y, P, (*Global::obs)(i,3), Elev, sigma_el, dEdY, 6);

        // Range and partials
        for (int k = 1; k <=3; k++){
            r(k) = Y(k);
        }

        s = transpose(LT*(U*transpose(r)-transpose(Rs)));                          // Topocentric position [m]

        Dist = norm(s);
        dDds = s/Dist;         // Range

        aux3 = dDds*LT*U;
        aux4 = zeros(1,3);

        for (int k = 1; k <= 3; k++){
            dDdY(k) = aux3(k);
            dDdY(k+3) = aux4(k);
        }

        // Measurement update
        MeasUpdate(K, Y, P, (*Global::obs)(i,4), Dist, sigma_range, dDdY, 6);
    }
    IERS(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, eop, (*Global::obs)(46,1), 'l');
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
    Global::auxparam.Mjd_UTC = Mjd_UTC;
    Global::auxparam.Mjd_TT = Mjd_TT;

    n_eqn = 6;
    Matrix Y0 = DEInteg (Accel, n_eqn, Y, 0, -((*Global::obs)(46,1)-(*Global::obs)(1,1))*86400.0, 1e-13, 1e-6);

    Matrix Y_true (6);
    Y_true(1) = 5753.173e3;     Y_true(2) = 2673.361e3;     Y_true(3) = 3440.304e3;
    Y_true(4) = 4.324207e3;     Y_true(5) = -1.924299e3;     Y_true(6) = -5.728216e3;

    printf("\nError of Position Estimation\n");
    printf("dX      %10.1f [m]\n",Y0(1)-Y_true(1));
    printf("dY      %10.1f [m]\n",Y0(2)-Y_true(2));
    printf("dZ      %10.1f [m]\n",Y0(3)-Y_true(3));
    printf("\nError of Velocity Estimation\n");
    printf("dVx     %8.1f [m/s]\n",Y0(4)-Y_true(4));
    printf("dVy     %8.1f [m/s]\n",Y0(5)-Y_true(5));
    printf("dVz     %8.1f [m/s]\n",Y0(6)-Y_true(6));

    return 0;
}
