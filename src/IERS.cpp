// $Source$
//----------------------------------------------------------------------------------------
//                          IERS
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/24
//
/*
 * @file IERS.cpp
 * @brief Management of IERS time and polar motion data
 *
 * @details This file contains declarations related to the management of IERS time and polar motion data.
 * @param x_pole Output parameter for x-coordinate of the Celestial Intermediate Pole (CIP) [rad].
 * @param y_pole Output parameter for y-coordinate of the Celestial Intermediate Pole (CIP) [rad].
 * @param UT1_UTC Output parameter for the difference between Universal Time (UT1) and Coordinated Universal Time (UTC) [s].
 * @param LOD Output parameter for the length of day [s].
 * @param dpsi Output parameter for the nutation correction in longitude [rad].
 * @param deps Output parameter for the nutation correction in obliquity [rad].
 * @param dx_pole Output parameter for the rate of change of the x-coordinate of the Celestial Intermediate Pole (CIP) [rad/s].
 * @param dy_pole Output parameter for the rate of change of the y-coordinate of the Celestial Intermediate Pole (CIP) [rad/s].
 * @param TAI_UTC Output parameter for the difference between International Atomic Time (TAI) and Coordinated Universal Time (UTC) [s].
 * @param eop Matrix containing Earth Orientation Parameters (EOP) data.
 * @param Mjd_UTC Modified Julian Date (UTC).
 * @param interp Interpolation type ('l' for linear, 'n' for nearest neighbor).
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\IERS.h"

void IERS(double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi,
          double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC, Matrix eop, 
          double Mjd_UTC, char interp)
{

    if (interp == 'l')
    {
        // linear interpolation
        double mjd = (floor(Mjd_UTC));
        int i;
        for(i = 1; i <= eop.n_column; i++){
            if(eop(4, i)==mjd){
                break;
            }
        }
        Matrix preeop = transpose(extract_column(eop, i));
        Matrix nexteop = transpose(extract_column(eop, i+1));
        double mfme = 1440.0 * (Mjd_UTC - floor(Mjd_UTC));
        double fixf = mfme / 1440.0;
        // Setting of IERS Earth rotation parameters 
        // (UT1 - UTC[s], TAI - UTC[s], x["], y ["]) 
        x_pole = preeop(5) + (nexteop(5) - preeop(5)) * fixf;
        y_pole = preeop(6) + (nexteop(6) - preeop(6)) * fixf;
        UT1_UTC = preeop(7) + (nexteop(7) - preeop(7)) * fixf;
        LOD = preeop(8) + (nexteop(8) - preeop(8)) * fixf;
        dpsi = preeop(9) + (nexteop(9) - preeop(9)) * fixf;
        deps = preeop(10) + (nexteop(10) - preeop(10)) * fixf;
        dx_pole = preeop(11) + (nexteop(11) - preeop(11)) * fixf;
        dy_pole = preeop(12) + (nexteop(12) - preeop(12)) * fixf;
        TAI_UTC = preeop(13);
        
        x_pole = x_pole / Const::Arcs; // Pole coordinate[rad] 
        y_pole = y_pole / Const::Arcs; // Pole coordinate[rad] 
        dpsi = dpsi / Const::Arcs;
        deps = deps / Const::Arcs;
        dx_pole = dx_pole / Const::Arcs; // Pole coordinate[rad] 
        dy_pole = dy_pole / Const::Arcs; // Pole coordinate[rad]
    }
    else
    {
        if (interp == 'n')
        {
            double mjd = (floor(Mjd_UTC));
            int i;
            for(i = 1; i <= eop.n_column; i++){
                if(eop(4, i)==mjd){
                    break;
                }
            }

            eop = transpose(extract_column(eop, i));
            // Setting of IERS Earth rotation parameters 
            // (UT1 - UTC[s], TAI - UTC[s], x["], y ["]) 
            x_pole = eop(5) / Const::Arcs; // Pole coordinate[rad] 
            y_pole = eop(6) / Const::Arcs; // Pole coordinate[rad] 
            UT1_UTC = eop(7); // UT1 - UTC time difference[s] 
            LOD = eop(8); // Length of day[s] 
            dpsi = eop(9) / Const::Arcs;
            deps = eop(10) / Const::Arcs;
            dx_pole = eop(11) / Const::Arcs; // Pole coordinate[rad] 
            dy_pole = eop(12) / Const::Arcs; // Pole coordinate[rad] 
            TAI_UTC = eop(13); // TAI - UTC time difference[s]
        }
    }
}