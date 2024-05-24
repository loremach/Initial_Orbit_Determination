//
// Created by lore0 on 23/05/2024.
//

#include "..\include\anglesg.h"

void anglesg (double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix Rs1, Matrix Rs2, Matrix Rs3, Matrix & r2, Matrix & v2){
    Matrix L1(3);   L1(1) = cos(el1)*sin(az1);  L1(2) = cos(el1)*cos(az1);  L1(3) = sin(el1);
    Matrix L2(3);   L2(1) = cos(el2)*sin(az2);  L2(2) = cos(el2)*cos(az2);  L2(3) = sin(el2);
    Matrix L3(3);   L3(1) = cos(el3)*sin(az3);  L3(2) = cos(el3)*cos(az3);  L3(3) = sin(el3);

    double lon1, lat1, h1, lon2, lat2, h2, lon3, lat3, h3;
    Geodetic(Rs1, lon1, lat1, h1);
    Geodetic(Rs2, lon2, lat2, h2);
    Geodetic(Rs3, lon3, lat3, h3);

    Matrix M1 = LTC(lon1, lat1);
    Matrix M2 = LTC(lon2, lat2);
    Matrix M3 = LTC(lon3, lat3);

    // body-fixed system
    Matrix Lb1 = transpose(M1)*transpose(L1);
    Matrix Lb2 = transpose(M1)*transpose(L2);
    Matrix Lb3 = transpose(M1)*transpose(L3);

    // mean of date system (J2000)
    double Mjd_UTC = Mjd1;
    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;
    IERS(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, *Global::eopdata, Mjd_UTC, 'l');

    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    double Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
    double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
    Matrix P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;

    Matrix pole = PoleMatrix(x_pole,y_pole);

    Matrix gha = GHAMatrix(Mjd_UT1);
    Matrix E = pole * gha * T;

    Matrix Lm1 = transpose(E)*Lb1;
    Rs1 = transpose(E)*transpose(Rs1);

    Mjd_UTC = Mjd2;
    IERS(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, *Global::eopdata, Mjd_UTC, 'l');
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    pole = PoleMatrix(x_pole,y_pole);
    gha = GHAMatrix(Mjd_UT1);
    E = pole * gha * T;

    Matrix Lm2 = transpose(E)*Lb2;
    Rs2 = transpose(E)* transpose(Rs2);

    Mjd_UTC = Mjd3;
    IERS(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, *Global::eopdata, Mjd_UTC, 'l');
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    pole = PoleMatrix(x_pole,y_pole);
    gha = GHAMatrix(Mjd_UT1);
    E = pole * gha * T;

    Matrix Lm3 = transpose(E)*Lb3;
    Rs3 = transpose(E)*transpose(Rs3);

    // geocentric inertial position
    double tau1 = (Mjd1-Mjd2)*86400.0;
    double tau3 = (Mjd3-Mjd2)*86400.0;

    double a1 = tau3/(tau3-tau1);
    double a3 =-tau1/(tau3-tau1);

    double b1 = tau3/(6*(tau3-tau1))*(pow(tau3-tau1,2)-pow(tau3,2));
    double b3 =-tau1/(6*(tau3-tau1))*(pow(tau3-tau1,2)-pow(tau1,2));

    Matrix aux1 = zeros(3,3);
    aux1 = assign_column(aux1, 1, Lm1);
    aux1 = assign_column(aux1, 2, Lm2);
    aux1 = assign_column(aux1, 3, Lm3);
    Matrix aux2 = zeros(3,3);
    aux2 = assign_column(aux2, 1, Rs1);
    aux2 = assign_column(aux2, 2, Rs2);
    aux2 = assign_column(aux2, 3, Rs2);

    Matrix D = inv(aux1)*aux2;

    double d1s = D(2,1)*a1-D(2,2)+D(2,3)*a3;
    double d2s = D(2,1)*b1+D(2,3)*b3;

    double Ccye = 2*dot(transpose(Lm2),transpose(Rs2));

    Matrix poly(9);
    poly(1)=  1.0;  // R2^8... polynomial
    poly(2)=  0.0;
    poly(3)=  -(pow(d1s,2) + d1s*Ccye + pow(norm(Rs2),2));
    poly(4)=  0.0;
    poly(5)=  0.0;
    poly(6)=  -Const::GM_Earth*(d2s*Ccye + 2*d1s*d2s);
    poly(7)=  0.0;
    poly(8)=  0.0;
    poly(9)=  -pow(Const::GM_Earth,2)*pow(d2s,2);

    Matrix rootsr(8);
    Matrix rootsi(8);
    roots(poly, rootsr, rootsi);

    double bigr2= -99999990.0;

    for (int j=1; j <=8; j++) {
        if (rootsr(j) > bigr2 && rootsi(j)==0){
            bigr2 = rootsr(j);
        }
    }
    double u = Const::GM_Earth/(pow(bigr2,3));

    double C1 = a1+b1*u;
    double C2 = -1.0;
    double C3 = a3+b3*u;

    Matrix auxC(3);     auxC(1) = C1;   auxC(2) = C2;   auxC(3) = C3;
    Matrix temp = D*transpose(auxC)*(-1);
    double rho1 = temp(1, 1)/(a1+b1*u);
    double rho2 = -temp(2, 1);
    double rho3 = temp(3, 1)/(a3+b3*u);

    double rhoold1 = rho1;
    double rhoold2 = rho2;
    double rhoold3 = rho3;

    rho2 = 99999999.9;
    double ll   = 0;

    while ((fabs(rhoold2-rho2) > 1e-12) && (ll <= 0 )) {
        ll = ll + 1;
        rho2 = rhoold2;

        Matrix r1 = Rs1 + Lm1 * rho1;
        r2 = Rs2 + Lm2 * rho2;
        Matrix r3 = Rs3 + Lm3* rho3;

        double magr1 = norm(r1);
        double magr2 = norm(r2);
        double magr3 = norm(r3);

        double theta, theta1, cop;
        string error = "";

        gibbs(transpose(r1), transpose(r2), transpose(r3), v2, theta, theta1, cop, error);

        if ((error.compare("ok") != 0) & (cop < Const::pi / 180.0)) {
            hgibbs(transpose(r1), transpose(r2), transpose(r3), Mjd1, Mjd2, Mjd3, v2, theta, theta1, cop, error);
        }

        Matrix aux3(6);
        for(int i = 1; i<=3; i++){
            aux3(i) = r2(i, 1);
            aux3(i+3) = v2(i);
        }
        double p;
        double a;
        double e;
        double i;
        double Omega;
        double omega;
        double M;
        elements(p, a, e, i, Omega, omega, M, aux3);
        double f1, g1, f3, g3;
        if (ll <= 8) {
            u = Const::GM_Earth / pow(magr2,3);

            double rdot = dot(transpose(r2), v2) / magr2;
            double udot = (-3 * Const::GM_Earth * rdot) / (pow(magr2, 4));

            double tausqr = tau1 * tau1;

            f1 = 1.0 - 0.5 * u * tausqr - (1.0 / 6.0) * udot * tausqr * tau1
                        - (1.0 / 24.0) * u * u * tausqr * tausqr
                        - (1.0 / 30.0) * u * udot * tausqr * tausqr * tau1;
            g1 = tau1 - (1.0 / 6.0) * u * tau1 * tausqr - (1.0 / 12.0) * udot * tausqr * tausqr
                        - (1.0 / 120.0) * u * u * tausqr * tausqr * tau1
                        - (1.0 / 120.0) * u * udot * tausqr * tausqr * tausqr;
            tausqr = tau3 * tau3;
            f3 = 1.0 - 0.5 * u * tausqr - (1.0 / 6.0) * udot * tausqr * tau3
                        - (1.0 / 24.0) * u * u * tausqr * tausqr
                        - (1.0 / 30.0) * u * udot * tausqr * tausqr * tau3;
            g3 = tau3 - (1.0 / 6.0) * u * tau3 * tausqr - (1.0 / 12.0) * udot * tausqr * tausqr
                        - (1.0 / 120.0) * u * u * tausqr * tausqr * tau3
                        - (1.0 / 120.0) * u * udot * tausqr * tausqr * tausqr;

        }else {

            theta = angl(transpose(r1), transpose(r2));
            theta1 = angl(transpose(r2), transpose(r3));

            f1 = 1 - ((magr1 * (1 - cos(theta)) / p));
            g1 = (magr1 * magr2 * sin(-theta)) / sqrt(p);
            f3 = 1 - ((magr3 * (1 - cos(theta1)) / p));
            g3 = (magr3 * magr2 * sin(theta1)) / sqrt(p);
        }

        C1 = g3 / (f1 * g3 - f3 * g1);
        C2 = -1;
        C3 = -g1 / (f1 * g3 - f3 * g1);

        double H1 = Const::GM_Earth * tau3 / 12;
        double H3 = -Const::GM_Earth * tau1 / 12;
        double H2 = H1 - H3;

        double G1 = -tau3 / (tau1 * (tau3 - tau1));
        double G3 = -tau1 / (tau3 * (tau3 - tau1));
        double G2 = G1 - G3;

        double D1 = G1 + H1 / pow(magr1, 3);
        double D2 = G2 + H2 / pow(magr2, 3);
        double D3 = G3 + H3 / pow(magr3, 3);

        Matrix auxD(3);     auxD(1) = D1;   auxD(2) = D2;   auxD(3) = D3;
        Matrix auxC2(3);     auxC2(1) = C1;   auxC2(2) = C2;   auxC2(3) = C3;
        temp = auxD*transpose(auxC2)*(-1);

        rhoold1 = temp(1) / (a1 + b1 * u);
        rhoold2 = -temp(1);
        rhoold3 = temp(1) / (a3 + b3 * u);

        r1 = Rs1 + Lm1 * rhoold1;
        r2 = Rs2 + Lm2 * rhoold2;
        r3 = Rs3 + Lm3 * rhoold3;
    }
    Matrix r1 = transpose(Rs1+Lm1*rho1);
    r2 = transpose(Rs2+Lm2*rho2);
    Matrix r3 = transpose(Rs3+Lm3*rho3);

}