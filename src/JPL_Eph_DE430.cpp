// $Source$
//----------------------------------------------------------------------------------------
//                          JPL_Eph_DE430
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/10
//
/*
 * @file JPL_Eph_DE430.cpp
 * @brief Computes the positions of the sun, moon, and nine major planets' equatorial positions using JPL Ephemerides
 *
 * @details This file contains the implementation of the function to compute the positions of the sun, moon, and nine major planets' equatorial positions using JPL Ephemerides.
 *
 * @param Mjd_TDB Modified Julian Date of TDB
 * @param[out] r_Earth Solar system barycenter (SSB) position
 * @param[out] r_Mars Mars position
 * @param[out] r_Mercury Mercury position
 * @param[out] r_Venus Venus position
 * @param[out] r_Jupiter Jupiter position
 * @param[out] r_Saturn Saturn position
 * @param[out] r_Uranus Uranus position
 * @param[out] r_Neptune Neptune position
 * @param[out] r_Pluto Pluto position
 * @param[out] r_Moon Moon position
 * @param[out] r_Sun Geocentric equatorial position referred to the International Celestial Reference Frame (ICRF)
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\JPL_Eph_DE430.h"

void JPL_Eph_DE430(Matrix &r_Mercury, Matrix &r_Venus, Matrix &r_Earth, Matrix &r_Mars,
                   Matrix &r_Jupiter, Matrix &r_Saturn, Matrix &r_Uranus, Matrix &r_Neptune,
                   Matrix &r_Pluto, Matrix &r_Moon, Matrix &r_Sun, double Mjd_TDB)
{
    double JD = Mjd_TDB + 2400000.5;
    int m;
    for (m = 1; m <= 2285; m++)
    {
        if ((*Global::PC)(m, 1) <= JD && JD <= (*Global::PC)(m, 2))
        {
            break;
        }
    }
    
    // i = find(PC(:,1)<=JD & JD<=PC(:,2),1,'first');
    Matrix PCtemp = extract_row((*Global::PC), m);
    double t1 = PCtemp(1) - 2400000.5; // MJD at start of interval
    double dt = Mjd_TDB - t1;
    Matrix temp(4);
    int k = 1;
    for (int l = 231; l <= 270; l = l + 13)
    {
        temp(k) = l;
        k++;
    }
    // temp = (231:13:270); //Empieza en 232, termina en 270 y salta de 13 en 13
    Matrix Cx_Earth(26);
    Matrix Cy_Earth(26);
    Matrix Cz_Earth(26);
    for (int n = 1; n <= 13; n++)
    {
        Cx_Earth(n) = PCtemp(temp(1) + n - 1);
        Cy_Earth(n) = PCtemp(temp(2) + n - 1);
        Cz_Earth(n) = PCtemp(temp(3) + n - 1);
    }
    
    temp = temp + 39;

    Matrix Cx(13);
    Matrix Cy(13);
    Matrix Cz(13);

    for (int n = 1; n <= 13; n++)
    {
        Cx(n) = PCtemp(temp(1) + n - 1);
        Cy(n) = PCtemp(temp(2) + n - 1);
        Cz(n) = PCtemp(temp(3) + n - 1);
    }

    for (int n = 1; n <= 13; n++)
    {
        Cx_Earth(n + 13) = Cx(n);
        Cy_Earth(n + 13) = Cy(n);
        Cz_Earth(n + 13) = Cz(n);
    }

    int j = -1;
    double Mjd0 = 0.0;
    if (0 <= dt && dt <= 16)
    {
        j = 0;
        Mjd0 = t1;
    }
    else if (16 < dt && dt <= 32)
    {
        j = 1;
        Mjd0 = t1 + 16.0 * j;
    }

    Matrix cheb_x(13);
    Matrix cheb_y(13);
    Matrix cheb_z(13);

    for (int n = 1; n <= 13; n++)
    {
        cheb_x(n) = Cx_Earth(13 * j + n);
        cheb_y(n) = Cy_Earth(13 * j + n);
        cheb_z(n) = Cz_Earth(13 * j + n);
    }

    Matrix Cheb = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 16, cheb_x, cheb_y, cheb_z);
    r_Earth = Cheb * 1e3;

    k = 1;
    for (int l = 441; l <= 480; l = l + 13)
    {
        temp(k) = l;
        k++;
    }

    Matrix Cx_Moon(109);
    Matrix Cy_Moon(109);
    Matrix Cz_Moon(109);
    for (int n = 1; n <= 13; n++)
    {
        Cx_Moon(n) = PCtemp(temp(1) + n - 1);
        Cy_Moon(n) = PCtemp(temp(2) + n - 1);
        Cz_Moon(n) = PCtemp(temp(3) + n - 1);
    }

    for (int i = 1; i <= 7; i++)
    {
        temp = temp + 39;
        for (int n = 1; n <= 13; n++)
        {
            Cx(n) = PCtemp(temp(1) + n - 1);
            Cy(n) = PCtemp(temp(2) + n - 1);
            Cz(n) = PCtemp(temp(3) + n - 1);
        }

        for (int n = 1; n <= 13; n++)
        {
            Cx_Moon(n + 13*i) = Cx(n);
            Cy_Moon(n + 13*i) = Cy(n);
            Cz_Moon(n + 13*i) = Cz(n);
        }
    }

    if (0 <= dt && dt <= 4){
        j = 0;
        Mjd0 = t1;
    }
    else if(4 < dt && dt <= 8){
        j = 1;
        Mjd0 = t1 + 4 * j;
    }else if(8 < dt && dt <= 12){
        j = 2;
        Mjd0 = t1 + 4 * j;
    }else if(12 < dt && dt <= 16){
        j = 3;
        Mjd0 = t1 + 4 * j;
    }else if(16 < dt && dt <= 20){
        j = 4;
        Mjd0 = t1 + 4 * j;
    }else if(20 < dt && dt <= 24){
        j = 5;
        Mjd0 = t1 + 4 * j;
    }else if(24 < dt && dt <= 28){
        j = 6;
        Mjd0 = t1 + 4 * j;
    }else if(28 < dt && dt <= 32){
        j = 7;
        Mjd0 = t1 + 4 * j;
    }
    
    for (int n = 1; n <= 13; n++)
    {
        cheb_x(n) = Cx_Moon(13 * j + n);
        cheb_y(n) = Cy_Moon(13 * j + n);
        cheb_z(n) = Cz_Moon(13 * j + n);
    }
    Cheb = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 4, cheb_x, cheb_y, cheb_z);
    r_Moon = Cheb*1e3;

    k = 1;
    for (int l = 753; l <= 786; l = l + 11)
    {
        temp(k) = l;
        k++;
    }

    Matrix Cx_Sun(22);
    Matrix Cy_Sun(22);
    Matrix Cz_Sun(22);

    for (int n = 1; n <= 11; n++)
    {
        Cx_Sun(n) = PCtemp(temp(1) + n - 1);
        Cy_Sun(n) = PCtemp(temp(2) + n - 1);
        Cz_Sun(n) = PCtemp(temp(3) + n - 1);
    }

    temp = temp + 33;
    Cx = zeros(1, 11);
    Cy = zeros(1, 11);
    Cz = zeros(1, 11);
    // Matrix Cx(11);
    // Matrix Cy(11);
    // Matrix Cz(11);

    for (int n = 1; n <= 11; n++)
    {
        Cx(n) = PCtemp(temp(1) + n - 1);
        Cy(n) = PCtemp(temp(2) + n - 1);
        Cz(n) = PCtemp(temp(3) + n - 1);
    }

    for (int n = 1; n <= 11; n++)
    {
        Cx_Sun(n + 11) = Cx(n);
        Cy_Sun(n + 11) = Cy(n);
        Cz_Sun(n + 11) = Cz(n);
    }
    
    if (0 <= dt && dt <= 16){
        j = 0;
        Mjd0 = t1;
    }else if(16 < dt && dt <= 32){
        j = 1;
        Mjd0 = t1 + 16 * j;
    }
    
    cheb_x = zeros(1, 11);
    cheb_y = zeros(1, 11);
    cheb_z = zeros(1, 11);
    // Matrix cheb_x(11);
    // Matrix cheb_y(11);
    // Matrix cheb_z(11);
    
    for (int n = 1; n <= 11; n++)
    {
        cheb_x(n) = Cx_Sun(11 * j + n);
        cheb_y(n) = Cy_Sun(11 * j + n);
        cheb_z(n) = Cz_Sun(11 * j + n);
    }

    Cheb = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 16, cheb_x, cheb_y, cheb_z);
    r_Sun = Cheb * 1e3;
    
    k = 1;
    for (int l = 3; l <= 45; l = l + 14)
    {
        temp(k) = l;
        k++;
    }
    // temp = (3 : 14 : 45);

    Matrix Cx_Mercury(56);
    Matrix Cy_Mercury(56);
    Matrix Cz_Mercury(56);

    for (int n = 1; n <= 14; n++)
    {
        Cx_Mercury(n) = PCtemp(temp(1) + n - 1);
        Cy_Mercury(n) = PCtemp(temp(2) + n - 1);
        Cz_Mercury(n) = PCtemp(temp(3) + n - 1);
    }

    Cx = zeros(1, 14);
    Cy = zeros(1, 14);
    Cz = zeros(1, 14);

    // Matrix Cx(14);
    // Matrix Cy(14);
    // Matrix Cz(14);

    for (int i = 1; i <= 3; i++)
    {
        temp = temp + 42;
        for (int n = 1; n <= 14; n++)
        {
            Cx(n) = PCtemp(temp(1) + n - 1);
            Cy(n) = PCtemp(temp(2) + n - 1);
            Cz(n) = PCtemp(temp(3) + n - 1);
        }

        for (int n = 1; n <= 14; n++)
        {
            Cx_Mercury(n + 14*i) = Cx(n);
            Cy_Mercury(n + 14*i) = Cy(n);
            Cz_Mercury(n + 14*i) = Cz(n);
        }
    }

    if (0 <= dt && dt <= 8){
        j = 0;
        Mjd0 = t1;
    }else if(8 < dt && dt <= 16){
        j = 1;
        Mjd0 = t1 + 8 * j;
    }else if(16 < dt && dt <= 24){
        j = 2;
        Mjd0 = t1 + 8 * j;
    }else if(24 < dt && dt <= 32){
        j = 3;
        Mjd0 = t1 + 8 * j;
    }

    cheb_x = zeros(1, 14);
    cheb_y = zeros(1, 14);
    cheb_z = zeros(1, 14);

    // Matrix cheb_x(14);
    // Matrix cheb_y(14);
    // Matrix cheb_z(14);
    
    for (int n = 1; n <= 14; n++)
    {
        cheb_x(n) = Cx_Mercury(14 * j + n);
        cheb_y(n) = Cy_Mercury(14 * j + n);
        cheb_z(n) = Cz_Mercury(14 * j + n);
    }

    Cheb = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0 + 8, cheb_x, cheb_y, cheb_z);
    r_Mercury = Cheb * 1e3;
    
    k = 1;
    for (int l = 171; l <= 201; l = l + 10)
    {
        temp(k) = l;
        k++;
    }

    // temp = (171 : 10 : 201);

    Matrix Cx_Venus(20);
    Matrix Cy_Venus(20);
    Matrix Cz_Venus(20);

    for (int n = 1; n <= 10; n++)
    {
        Cx_Venus(n) = PCtemp(temp(1) + n - 1);
        Cy_Venus(n) = PCtemp(temp(2) + n - 1);
        Cz_Venus(n) = PCtemp(temp(3) + n - 1);
    }

    temp = temp + 30;

    Cx = zeros(1, 10);
    Cy = zeros(1, 10);
    Cz = zeros(1, 10);

    // Matrix Cx(10);
    // Matrix Cy(10);
    // Matrix Cz(10);

    for (int n = 1; n <= 10; n++)
    {
        Cx(n) = PCtemp(temp(1) + n - 1);
        Cy(n) = PCtemp(temp(2) + n - 1);
        Cz(n) = PCtemp(temp(3) + n - 1);
    }

    for (int n = 1; n <= 10; n++)
    {
        Cx_Venus(n + 10) = Cx(n);
        Cy_Venus(n + 10) = Cy(n);
        Cz_Venus(n + 10) = Cz(n);
    }

    if (0 <= dt && dt <= 16){
        j = 0;
        Mjd0 = t1;
    }else if(16 < dt && dt <= 32){
        j = 1;
        Mjd0 = t1 + 16 * j;
    }

    cheb_x = zeros(1, 10);
    cheb_y = zeros(1, 10);
    cheb_z = zeros(1, 10);
    // Matrix cheb_x(10);
    // Matrix cheb_y(10);
    // Matrix cheb_z(10);
    
    for (int n = 1; n <= 10; n++)
    {
        cheb_x(n) = Cx_Venus(10 * j + n);
        cheb_y(n) = Cy_Venus(10 * j + n);
        cheb_z(n) = Cz_Venus(10 * j + n);
    }

    Cheb = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 16, cheb_x, cheb_y, cheb_z);
    r_Venus = Cheb * 1e3;
    
    k = 1;
    for (int l = 309; l <= 342; l = l + 11)
    {
        temp(k) = l;
        k++;
    }
    // temp = (309 : 11 : 342);
    Matrix Cx_Mars(22);
    Matrix Cy_Mars(22);
    Matrix Cz_Mars(22);

    for (int n = 1; n <= 11; n++)
    {
        Cx_Mars(n) = PCtemp(temp(1) + n - 1);
        Cy_Mars(n) = PCtemp(temp(2) + n - 1);
        Cz_Mars(n) = PCtemp(temp(3) + n - 1);
    }

    j = 0;
    Mjd0 = t1;
    cheb_x = zeros(1, 11);
    cheb_y = zeros(1, 11);
    cheb_z = zeros(1, 11);

    // Matrix cheb_x(11);
    // Matrix cheb_y(11);
    // Matrix cheb_z(11);
    
    for (int n = 1; n <= 11; n++)
    {
        cheb_x(n) = Cx_Mars(11 * j + n);
        cheb_y(n) = Cy_Mars(11 * j + n);
        cheb_z(n) = Cz_Mars(11 * j + n);
    }

    Cheb = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 32, cheb_x, cheb_y, cheb_z);
    r_Mars = Cheb * 1e3;
    
    k = 1;
    for (int l = 342; l <= 366; l = l + 8)
    {
        temp(k) = l;
        k++;
    }
    // temp = (342 : 8 : 366);

    Matrix Cx_Jupiter(16);
    Matrix Cy_Jupiter(16);
    Matrix Cz_Jupiter(16);

    for (int n = 1; n <= 8; n++)
    {
        Cx_Jupiter(n) = PCtemp(temp(1) + n - 1);
        Cy_Jupiter(n) = PCtemp(temp(2) + n - 1);
        Cz_Jupiter(n) = PCtemp(temp(3) + n - 1);
    }

    j = 0;
    Mjd0 = t1;

    cheb_x = zeros(1, 8);
    cheb_y = zeros(1, 8);
    cheb_z = zeros(1, 8);

    // Matrix cheb_x(8);
    // Matrix cheb_y(8);
    // Matrix cheb_z(8);
    
    for (int n = 1; n <= 8; n++)
    {
        cheb_x(n) = Cx_Jupiter(8 * j + n);
        cheb_y(n) = Cy_Jupiter(8 * j + n);
        cheb_z(n) = Cz_Jupiter(8 * j + n);
    }

    Cheb = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0 + 32, cheb_x, cheb_y, cheb_z);
    r_Jupiter = Cheb * 1e3;
    
    k = 1;
    for (int l = 366; l <= 387; l = l + 7)
    {
        temp(k) = l;
        k++;
    }
    // temp = (366 : 7 : 387);
    Matrix Cx_Saturn(14);
    Matrix Cy_Saturn(14);
    Matrix Cz_Saturn(14);

    for (int n = 1; n <= 7; n++)
    {
        Cx_Saturn(n) = PCtemp(temp(1) + n - 1);
        Cy_Saturn(n) = PCtemp(temp(2) + n - 1);
        Cz_Saturn(n) = PCtemp(temp(3) + n - 1);
    }

    j = 0;
    Mjd0 = t1;

    cheb_x = zeros(1, 7);
    cheb_y = zeros(1, 7);
    cheb_z = zeros(1, 7);

    // Matrix cheb_x(7);
    // Matrix cheb_y(7);
    // Matrix cheb_z(7);
    
    for (int n = 1; n <= 7; n++)
    {
        cheb_x(n) = Cx_Saturn(7 * j + n);
        cheb_y(n) = Cy_Saturn(7 * j + n);
        cheb_z(n) = Cz_Saturn(7 * j + n);
    }

    Cheb = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0 + 32, cheb_x, cheb_y, cheb_z);
    r_Saturn = Cheb * 1e3;

    k = 1;
    for (int l = 387; l <= 405; l = l + 6)
    {
        temp(k) = l;
        k++;
    }

    // temp = (387 : 6 : 405);

    Matrix Cx_Uranus(12);
    Matrix Cy_Uranus(12);
    Matrix Cz_Uranus(12);

    for (int n = 1; n <= 6; n++)
    {
        Cx_Uranus(n) = PCtemp(temp(1) + n - 1);
        Cy_Uranus(n) = PCtemp(temp(2) + n - 1);
        Cz_Uranus(n) = PCtemp(temp(3) + n - 1);
    }

    j = 0;
    Mjd0 = t1;

    cheb_x = zeros(1, 6);
    cheb_y = zeros(1, 6);
    cheb_z = zeros(1, 6);

    // Matrix cheb_x(6);
    // Matrix cheb_y(6);
    // Matrix cheb_z(6);
    
    for (int n = 1; n <= 6; n++)
    {
        cheb_x(n) = Cx_Uranus(6 * j + n);
        cheb_y(n) = Cy_Uranus(6 * j + n);
        cheb_z(n) = Cz_Uranus(6 * j + n);
    }

    Cheb = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, cheb_x, cheb_y, cheb_z);
    r_Uranus = Cheb * 1e3;
    
    k = 1;
    for (int l = 405; l <= 423; l = l + 6)
    {
        temp(k) = l;
        k++;
    }

    // temp = (405 : 6 : 423);

    Matrix Cx_Neptune(12);
    Matrix Cy_Neptune(12);
    Matrix Cz_Neptune(12);

    for (int n = 1; n <= 6; n++)
    {
        Cx_Neptune(n) = PCtemp(temp(1) + n - 1);
        Cy_Neptune(n) = PCtemp(temp(2) + n - 1);
        Cz_Neptune(n) = PCtemp(temp(3) + n - 1);
    }

    j = 0;
    Mjd0 = t1;
    
    cheb_x = zeros(1, 6);
    cheb_y = zeros(1, 6);
    cheb_z = zeros(1, 6);

    // Matrix cheb_x(6);
    // Matrix cheb_y(6);
    // Matrix cheb_z(6);
    
    for (int n = 1; n <= 6; n++)
    {
        cheb_x(n) = Cx_Neptune(6 * j + n);
        cheb_y(n) = Cy_Neptune(6 * j + n);
        cheb_z(n) = Cz_Neptune(6 * j + n);
    }

    Cheb = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, cheb_x, cheb_y, cheb_z);
    r_Neptune = Cheb * 1e3;
    
    k = 1;
    for (int l = 423; l <= 441; l = l + 6)
    {
        temp(k) = l;
        k++;
    }
    // temp = (423 : 6 : 441);

    Matrix Cx_Pluto(12);
    Matrix Cy_Pluto(12);
    Matrix Cz_Pluto(12);

    for (int n = 1; n <= 6; n++)
    {
        Cx_Pluto(n) = PCtemp(temp(1) + n - 1);
        Cy_Pluto(n) = PCtemp(temp(2) + n - 1);
        Cz_Pluto(n) = PCtemp(temp(3) + n - 1);
    }

    j = 0;
    Mjd0 = t1;
    
    cheb_x = zeros(1, 6);
    cheb_y = zeros(1, 6);
    cheb_z = zeros(1, 6);

    // Matrix cheb_x(6);
    // Matrix cheb_y(6);
    // Matrix cheb_z(6);
    
    for (int n = 1; n <= 6; n++)
    {
        cheb_x(n) = Cx_Pluto(6 * j + n);
        cheb_y(n) = Cy_Pluto(6 * j + n);
        cheb_z(n) = Cz_Pluto(6 * j + n);
    }

    Cheb = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, cheb_x, cheb_y, cheb_z);
    r_Pluto = Cheb * 1e3;
    
    k = 1;
    for (int l = 819; l <= 839; l = l + 10)
    {
        temp(k) = l;
        k++;
    }
    // temp = (819 : 10 : 839);
    // Nutations y Librations no se usan!!!!

    Matrix Cx_Nutations(40);
    Matrix Cy_Nutations(40);

    for (int n = 1; n <= 10; n++)
    {
        Cx_Nutations(n) = PCtemp(temp(1) + n - 1);
        Cy_Nutations(n) = PCtemp(temp(2) + n - 1);
    }
    Cx = zeros(1, 10);
    Cy = zeros(1, 10);
    // Matrix Cx(10);
    // Matrix Cy(10);

    for(int i = 1; i <= 3; i++) {
        temp = temp + 20;
        for (int n = 1; n <= 10; n++)
        {
            Cx(n) = PCtemp(temp(1) + n - 1);
            Cy(n) = PCtemp(temp(2) + n - 1);
        }
        for (int n = 1; n <= 10; n++)
        {
            Cx_Nutations(n + 10*i) = Cx(n);
            Cy_Nutations(n + 10*i) = Cy(n);
        }
    } 

    if (0 <= dt && dt <= 8){
        j = 0;
        Mjd0 = t1;
    }else if(8 < dt && dt <= 16){
        j = 1;
        Mjd0 = t1 + 8 * j;
    }else if(16 < dt && dt <= 24){
        j = 2;
        Mjd0 = t1 + 8 * j;
    }else if(24 < dt && dt <= 32){
        j = 3;
        Mjd0 = t1 + 8 * j;
    }

    cheb_x = zeros(1, 10);
    cheb_y = zeros(1, 10);

    // Matrix cheb_x(10);
    // Matrix cheb_y(10);
    
    for (int n = 1; n <= 10; n++)
    {
        cheb_x(n) = Cx_Nutations(10 * j + n);
        cheb_y(n) = Cy_Nutations(10 * j + n);
    }

    Cheb = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 8, cheb_x, cheb_y, zeros(1,10));
    Matrix Nutations = Cheb * 1e3;

    k = 1;
    for (int l = 899; l <= 929; l = l + 10)
    {
        temp(k) = l;
        k++;
    }

    // temp = (899 : 10 : 929);

    Matrix Cx_Librations(40);
    Matrix Cy_Librations(40);
    Matrix Cz_Librations(40);

    for (int n = 1; n <= 10; n++)
    {
        Cx_Librations(n) = PCtemp(temp(1) + n - 1);
        Cy_Librations(n) = PCtemp(temp(2) + n - 1);
        Cz_Librations(n) = PCtemp(temp(3) + n - 1);
    }
    Cx = zeros(1, 10);
    Cy = zeros(1, 10);
    Cz = zeros(1, 10);
    
    // Matrix Cx(10);
    // Matrix Cy(10);
    // Matrix Cz(10);

    for(int i = 1; i <= 3; i++){
        temp = temp + 30;
        for (int n = 1; n <= 10; n++)
        {
            Cx(n) = PCtemp(temp(1) + n - 1);
            Cy(n) = PCtemp(temp(2) + n - 1);
            Cz(n) = PCtemp(temp(3) + n - 1);
        }

        for (int n = 1; n <= 10; n++)
        {
            Cx_Librations(n + 10*i) = Cx(n);
            Cy_Librations(n + 10*i) = Cy(n);
            Cz_Librations(n + 10*i) = Cz(n);
        }
    }

    if (0 <= dt && dt <= 8){
        j = 0;
        Mjd0 = t1;
    }else if(8 < dt && dt <= 16){
        j = 1;
        Mjd0 = t1 + 8 * j;
    }else if(16 < dt && dt <= 24){
        j = 2;
        Mjd0 = t1 + 8 * j;
    }else if(24 < dt && dt <= 32){
        j = 3;
        Mjd0 = t1 + 8 * j;
    }
    cheb_x = zeros(1, 10);
    cheb_y = zeros(1, 10);
    cheb_z = zeros(1, 10);

    // Matrix cheb_x(10);
    // Matrix cheb_y(10);
    // Matrix cheb_z(10);
    
    for (int n = 1; n <= 10; n++)
    {
        cheb_x(n) = Cx_Librations(10 * j + n);
        cheb_y(n) = Cy_Librations(10 * j + n);
        cheb_z(n) = Cz_Librations(10 * j + n);
    }

    Cheb = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 8, cheb_x, cheb_y, cheb_z);
    Matrix Librations = Cheb * 1e3;

    double EMRAT = 81.30056907419062; // DE430
    double EMRAT1 = 1.0 / (1.0 + EMRAT);
    
    r_Earth = r_Earth - r_Moon * EMRAT1;
    r_Mercury = r_Mercury -r_Earth;
    r_Venus = r_Venus -r_Earth ;
    r_Mars =  r_Mars -r_Earth;
    r_Jupiter = r_Jupiter -r_Earth ;
    r_Saturn = r_Saturn -r_Earth;
    r_Uranus = r_Uranus -r_Earth;
    r_Neptune = r_Neptune -r_Earth;
    r_Pluto = r_Pluto -r_Earth;
    r_Sun = r_Sun -r_Earth;
}