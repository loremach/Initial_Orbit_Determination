#include "..\include\Accel.h"

Matrix Accel(double x, Matrix & Y){
    double x_pole; 
    double y_pole; 
    double UT1_UTC; 
    double LOD; 
    double dpsi;
    double deps; 
    double dx_pole; 
    double dy_pole; 
    double TAI_UTC;
    IERS(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, *Global::eopdata, Global::auxparam.Mjd_UTC + x/86400, 'l');

    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double Mjd_UT1 = Global::auxparam.Mjd_UTC + x/86400.0 + UT1_UTC/86400.0;
    double Mjd_TT = Global::auxparam.Mjd_UTC + x/86400.0 + TT_UTC/86400.0;
    Matrix P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;

    Matrix pole = PoleMatrix(x_pole,y_pole);
    Matrix gha = GHAMatrix(Mjd_UT1);

    Matrix E =  pole * gha * T;

    double MJD_TDB = Mjday_TDB(Mjd_TT);

    Matrix r_Mercury(3);
    Matrix r_Venus(3);
    Matrix r_Earth(3); 
    Matrix r_Mars(3);
    Matrix r_Jupiter(3); 
    Matrix r_Saturn(3);
    Matrix r_Uranus(3);
    Matrix r_Neptune(3);
    Matrix r_Pluto(3); 
    Matrix r_Moon(3); 
    Matrix r_Sun(3);

    JPL_Eph_DE430(r_Mercury, r_Venus, r_Earth, r_Mars,
                   r_Jupiter, r_Saturn, r_Uranus, r_Neptune,
                   r_Pluto, r_Moon, r_Sun, MJD_TDB);

    // Acceleration due to harmonic gravity field
    Matrix aux(3);
    for (int i = 1; i <= 3; i++){
        aux(i) = Y(i);
    }

    Matrix a = AccelHarmonic(aux, E, Global::auxparam.n, Global::auxparam.m);
    // Luni-solar perturbations

    if (Global::auxparam.sun){
        Matrix accel1 = AccelPointMass(aux,r_Sun,Const::GM_Sun);
        a = a + accel1;
    }

    if (Global::auxparam.moon){
        Matrix accel2 = AccelPointMass(aux,r_Moon,Const::GM_Moon);
        a = a + accel2;
    }
    // Planetary perturbations
    if (Global::auxparam.planets){
        Matrix accel3 = AccelPointMass(aux,r_Mercury,Const::GM_Mercury);
        a = a + accel3;
        accel3 = AccelPointMass(aux,r_Venus,Const::GM_Venus);
        a = a + accel3;
        accel3 = AccelPointMass(aux,r_Mars,Const::GM_Mars);
        a = a + accel3;
        accel3 = AccelPointMass(aux,r_Jupiter,Const::GM_Jupiter);
        a = a + accel3;
        accel3 = AccelPointMass(aux,r_Saturn,Const::GM_Saturn);
        a = a + accel3;
        accel3 = AccelPointMass(aux,r_Uranus,Const::GM_Uranus);
        a = a + accel3;
        accel3 = AccelPointMass(aux,r_Neptune,Const::GM_Neptune);  
        a = a + accel3;
        accel3 = AccelPointMass(aux,r_Pluto,Const::GM_Pluto);
        a = a + accel3;
    }
    Matrix dY(6);

    for(int i = 1; i<=3; i++){
        dY(i) = Y(i+3);
        dY(i+3) = a(i);
    }
    return dY;
}