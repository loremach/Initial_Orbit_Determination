#include "..\include\R_x.h"
#include "..\include\R_y.h"
#include "..\include\R_z.h"
#include "..\include\sign_.h"
#include "..\include\Frac.h"
#include "..\include\Legendre.h"
#include "..\include\timediff.h"
#include "..\include\unit.h"
#include "..\include\global.h"
#include "..\include\accelPointMass.h"
#include "..\include\AzElPa.h"
#include "..\include\Cheb3D.h"
#include "..\include\EccAnom.h"
#include "..\include\Geodetic.h"
#include "..\include\IERS.h"
#include "..\include\MeanObliquity.h"
#include "..\include\Mjday.h"
#include "..\include\Mjday_TDB.h"
#include "..\include\Position.h"
#include "..\include\NutAngles.h"
#include "..\include\EqnEquinox.h"
#include "..\include\NutMatrix.h"
#include "..\include\PoleMatrix.h"
#include "..\include\PrecMatrix.h"
#include "..\include\TimeUpdate.h"
#include "..\include\angl.h"
#include "..\include\gmst.h"
#include "..\include\LTC.h"
#include "..\include\MeasUpdate.h"
#include "..\include\elements.h"
#include "..\include\AccelHarmonic.h"
#include "..\include\gast.h"
#include "..\include\GHAMatrix.h"
#include "..\include\G_AccelHarmonic.h"
#include "..\include\JPL_Eph_DE430.h"
#include "..\include\Accel.h"
#include "..\include\gibbs.h"
#include "..\include\hgibbs.h"
#include "..\include\VarEqn.h"
#include "..\include\DEInteg.hpp"
#include "..\include\anglesg.h"
#include "..\include\Matrix.h"

#include <stdio.h>
#include <cmath>
#include <string>


int tests_run = 0;

#define FAIL() printf("\nFailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) {FAIL(); return 1;} } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(double **A, double **B, int f, int c, double p){
    int equal = 1;

    for(int i = 0; i < f; i++){
        for(int j = 0; j < c; j++){
            if(abs(A[i][j]-B[i][j]) > p){
                printf("%2.20lf %2.20lf\n", A[i][j], B[i][j]);
                equal = 0;
            }
        }
    }
    return equal;
}

int r_x_01(){

    Matrix aux = R_x(0.5);

    Matrix A(3,3);

    A(1,1) = 1;  A(1,2) = 0;        A(1,3) = 0;
    A(2,1) = 0;  A(2,2) = 0.8775;   A(2,3) = 0.4794;
    A(3,1) = 0;  A(3,2) = -0.4794;  A(3,3) = 0.8776;

    _assert(m_equals(A.data, aux.data, 3, 3, 1e-4));

    return 0;
}

int r_x_02(){

    Matrix aux = R_x(0.73);

    Matrix A(3,3);

    A(1,1) = 1;  A(1,2) =  0;       A(1,3) = 0;
    A(2,1) = 0;  A(2,2) = 0.7452;   A(2,3) = 0.6669;
    A(3,1) = 0;  A(3,2) = -0.6669;  A(3,3) = 0.7452;

    _assert(m_equals(A.data, aux.data, 3, 3, 1e-4));

    return 0;
}

int r_y_01(){

    Matrix aux = R_y(0.5);

    Matrix A(3,3);

    A(1,1) = 0.8776;    A(1,2) = 0;     A(1,3) = -0.4794;
    A(2,1) = 0;         A(2,2) = 1;     A(2,3) = 0;
    A(3,1) = 0.4794;    A(3,2) = 0;     A(3,3) = 0.8776;

    _assert(m_equals(A.data, aux.data, 3, 3, 1e-4));

    return 0;
}

int r_y_02(){

    Matrix aux = R_y(0.65);

    Matrix A(3,3);

    A(1,1) = 0.7961;    A(1,2) = 0;     A(1,3) = -0.6052;
    A(2,1) = 0;         A(2,2) = 1;     A(2,3) = 0;
    A(3,1) = 0.6052;    A(3,2) = 0;     A(3,3) = 0.7961;

    _assert(m_equals(A.data, aux.data, 3, 3, 1e-4));

    return 0;
}


int r_z_01(){

    Matrix aux = R_z(0.5);

    Matrix A(3,3);

    A(1,1) = 0.8776;    A(1,2) = 0.4794;     A(1,3) = 0;
    A(2,1) = -0.4794;   A(2,2) = 0.8776;     A(2,3) = 0;
    A(3,1) = 0;         A(3,2) = 0;          A(3,3) = 1;

    _assert(m_equals(A.data, aux.data, 3, 3, 1e-4));

    return 0;
}

int r_z_02(){

    Matrix aux = R_z(0.43);

    Matrix A(3,3);

    A(1,1) = 0.9090;    A(1,2) = 0.4169;     A(1,3) = 0;
    A(2,1) = -0.4169;   A(2,2) = 0.9090;     A(2,3) = 0;
    A(3,1) = 0;         A(3,2) = 0;          A(3,3) = 1;

    _assert(m_equals(A.data, aux.data, 3, 3, 1e-4));

    return 0;
}

int sign_01(){
    double a = - 3.5;
    double b = 4.43;
    double res = sign_(a,b);
    _assert(res==3.5);

    return 0;
}

int sign_02(){
    double a = - 3.5;
    double b = - 4.43;
    double res = sign_(a,b);
    _assert(res==-3.5);

    return 0;
}

int Frac_01(){
    double a = - 3.5;
    double res = Frac(a);
    _assert(res==-0.5);

    return 0;
}

int Frac_02(){
    double a = 4.4335;
    double res = Frac(a);
    _assert(res==0.4335);

    return 0;
}

int legendre_01(){
    Matrix aux1(4,4);
    Matrix aux2(4,4);

    Matrix A(4,4);
    Matrix B(4,4);

    A(1,1) = 1;         A(1,2) = 0;         A(1,3) = 0;         A(1,4) = 0;
    A(2,1) = 1.5024;    A(2,2) = 0.8618;    A(2,3) = 0;         A(2,4) = 0;
    A(3,1) = 1.4057;    A(3,2) = 1.6716;    A(3,3) = 0.4794;    A(3,4) = 0;
    A(4,1) = 0.8745;    A(4,2) = 2.2267;    A(4,3) = 1.1003;    A(4,4) = 0.2577;

    B(1,1) = 0;         B(1,2) = 0;         B(1,3) = 0;         B(1,4) = 0;
    B(2,1) = 0.8618;    B(2,2) = -1.5024;   B(2,3) = 0;         B(2,4) = 0;
    B(3,1) = 2.8953;    B(3,2) = -1.9553;   B(3,3) = -1.6716;   B(3,4) = 0;
    B(4,1) = 5.4543;    B(4,2) = -0.4024;   B(4,3) = -3.2051;   B(4,4) = -1.3476;

    Legendre(3, 3, 1.05, aux1, aux2);
    
    _assert(m_equals(A.data, aux1.data, 4, 4, 1e-4));
    _assert(m_equals(B.data, aux2.data, 4, 4, 1e-4));

    return 0;
}


int legendre_02(){
    Matrix aux1(3,3);
    Matrix aux2(3,3);

    Matrix A(3,3);
    Matrix B(3,3);

    A(1,1) = 1;         A(1,2) = 0;         A(1,3) = 0;         
    A(2,1) = 1.4385;    A(2,2) = 0.9648;    A(2,3) = 0;         
    A(3,1) = 1.1954;    A(3,2) = 1.7917;    A(3,3) = 0.6008;   

    B(1,1) = 0;         B(1,2) = 0;         B(1,3) = 0;        
    B(2,1) = 0.9648;    B(2,2) = -1.4385;   B(2,3) = 0;         
    B(3,1) = 3.1033;    B(3,2) = -1.4696;   B(3,3) = -1.7917;   

    Legendre(2, 1, 0.98, aux1, aux2);

    _assert(m_equals(A.data, aux1.data, 3, 2, 1e-4));
    _assert(m_equals(B.data, aux2.data, 3, 2, 1e-4));

    return 0;
}

int legendre_03(){
    Matrix aux1(3,4);
    Matrix aux2(3,4);

    Matrix A(3,4);
    Matrix B(3,4);

    A(1,1) = 1;         A(1,2) = 0;         A(1,3) = 0;         A(1,4) = 0;
    A(2,1) = 1.6821;    A(2,2) = 0.4131;    A(2,3) = 0;         A(2,4) = 0;
    A(3,1) = 2.0453;    A(3,2) = 0.8970;    A(3,3) = 0.1101;    A(3,4) = 0;

    B(1,1) = 0;         B(1,2) = 0;         B(1,3) = 0;         B(1,4) = 0;
    B(2,1) = 0.4131;    B(2,2) = -1.6821;   B(2,3) = 0;         B(2,4) = 0;
    B(3,1) = 1.5536;    B(3,2) = -3.4325;   B(3,3) = -0.8970;   B(3,4) = 0;


    Legendre(2, 3, 1.33, aux1, aux2);
    
    _assert(m_equals(A.data, aux1.data, 3, 4, 1e-4));
    _assert(m_equals(B.data, aux2.data, 3, 4, 1e-4));

    return 0;
}


int timediff_01(){
    double UT1_UTC = 12;
    double TAI_UTC = -4;
    
    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS; 
    double TT_UTC;
    double GPS_UTC;

    double UT1_TAI_m = 16;
    double UTC_GPS_m = 23;
    double UT1_GPS_m = 35; 
    double TT_UTC_m = 28.1840;
    double GPS_UTC_m = -23;

    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double eps = 1e-4;

    _assert(abs(UT1_TAI-UT1_TAI_m)<eps && abs(UTC_GPS-UTC_GPS_m)<eps && abs(UT1_GPS-UT1_GPS_m)<eps 
            && abs(TT_UTC-TT_UTC_m)<eps && abs(GPS_UTC-GPS_UTC_m)<eps);

    return 0;
}

int timediff_02(){
    double UT1_UTC = -4.456;
    double TAI_UTC = -21.5;
    
    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS; 
    double TT_UTC;
    double GPS_UTC;

    double UT1_TAI_m = 17.044;
    double UTC_GPS_m = 40.5;
    double UT1_GPS_m = 36.044; 
    double TT_UTC_m = 10.684;
    double GPS_UTC_m = -40.5;

    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double eps = 1e-4;

    _assert(abs(UT1_TAI-UT1_TAI_m)<eps && abs(UTC_GPS-UTC_GPS_m)<eps && abs(UT1_GPS-UT1_GPS_m)<eps 
            && abs(TT_UTC-TT_UTC_m)<eps && abs(GPS_UTC-GPS_UTC_m)<eps);

    return 0;
}


int unit_01(){

    Matrix vec(3);
    Matrix outvec(3);
    Matrix result(3);

    vec(1) = 1;             vec(2) = 2;             vec(3) = 3;
    result(1) = 0.2673;     result(2) = 0.5345;     result(3) = 0.8018;
    outvec = unit(vec);

    _assert(m_equals(result.data, outvec.data, 1, 3, 1e-4));

    return 0;
}

int unit_02(){
    Matrix vec(3);
    Matrix outvec(3);
    Matrix result(3);
    
    vec(1) = -2.3;           vec(2) = 5.65;          vec(3) = 7.86;
    result(1) = -0.2312;     result(2) = 0.5679;     result(3) = 0.7900;

    outvec = unit(vec);
    
    _assert(m_equals(result.data, outvec.data, 1, 3, 1e-4));

    return 0;
}


int accelPointMass_01(){
    Matrix result(3);
    Matrix r(3);        r(1) = 1;       r(2) = 2;       r(3) = 3;
    Matrix s(3);        s(1) = 4;       s(2) = 5;       s(3) = 6;

    Matrix compare(3);              compare(1) = 0.007731656678682;
    compare(2) = 0.006991652935438; compare(3) = 0.006251649192193;

    result = AccelPointMass(r, s, 0.5);
    _assert(m_equals(compare.data, result.data, 1, 3, 1e-14));

    return 0;
}

int accelPointMass_02(){
    Matrix result(3);
    Matrix r(3);        r(1) = 5.6;       r(2) = 2.4;       r(3) = -8.5;
    Matrix s(3);        s(1) = 12.1;      s(2) = -4;        s(3) = 6.98;

    Matrix compare(3);              compare(1) = -0.002456577081167;
    compare(2) = 0.000174484669332; compare(3) = 0.000342244099929;

    result = AccelPointMass(r, s, 0.87);
    _assert(m_equals(compare.data, result.data, 1, 3, 1e-14));

    return 0;
}

int AzElPa_01(){
    Matrix dAds(3);
    Matrix dEds(3);  
    double Az;  double El;    
    Matrix s(3);        s(1) = 1;      s(2) = 2;        s(3) = 3;

    Matrix dAds_comp(3);  dAds_comp(1) = 0.4;   dAds_comp(2) = -0.2;   dAds_comp(3) = 0.0;
    Matrix dEds_comp(3);  dEds_comp(1) = -0.095831484749991;   dEds_comp(2) = -0.191662969499982;   dEds_comp(3) = 0.159719141249985;
    double Az_comp = 0.463647609000806;  double El_comp = 0.930274014115472;

    AzElPa(s, Az, El, dAds, dEds);
    
    _assert(m_equals(dAds_comp.data, dAds.data, 1, 3, 1e-14));
    _assert(m_equals(dEds_comp.data, dEds.data, 1, 3, 1e-14));
    _assert(abs(Az-Az_comp)<1e-5);
    _assert(abs(El-El_comp)<1e-5);

    return 0;
}


int AzElPa_02(){
    Matrix dAds(3);
    Matrix dEds(3);  
    double Az;  double El;    
    Matrix s(3);        s(1) = -7.87;      s(2) = -3.21;        s(3) = 3.66;

    Matrix dAds_comp(3);  dAds_comp(1) = -0.044434600849933;   dAds_comp(2) = 0.108940906133636;    dAds_comp(3) = 0.0;
    Matrix dEds_comp(3);  dEds_comp(1) = 0.039573505886700;    dEds_comp(2) = 0.016141163138031;    dEds_comp(3) = 0.099250443989455;
    double Az_comp = 4.325109711424560;  double El_comp = 0.406617021214728;

    AzElPa(s, Az, El, dAds, dEds);

    _assert(m_equals(dAds_comp.data, dAds.data, 1, 3, 1e-14));
    _assert(m_equals(dEds_comp.data, dEds.data, 1, 3, 1e-14));
    _assert(abs(Az-Az_comp)<1e-5);
    _assert(abs(El-El_comp)<1e-5);

    return 0;
}

int Cheb3D_01(){

    double t = 10.3; 
    int N = 3; 
    double Ta = 0;
    double Tb = 12.3;
    Matrix Cx(3); Matrix Cy(3); Matrix Cz(3);
    Matrix result(3);
    Matrix compare(3);  
    
    Cx(1) = 1;      Cx(2) = 2;        Cx(3) = 3;
    Cy(1) = 4;      Cy(2) = 5;        Cy(3) = 6;
    Cz(1) = 7;      Cz(2) = 8;        Cz(3) = 9;

    compare(1) = 2.081697402339877;   compare(2) = 6.838191552647235;    compare(3) = 11.594685702954592;

    result = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);

    _assert(m_equals(compare.data, result.data, 1, 3, 1e-14));

    return 0;
}

int Cheb3D_02(){

    double t = 32.3; 
    int N = 3; 
    double Ta = 3.4;
    double Tb = 44.3;
    Matrix Cx(3); Matrix Cy(3); Matrix Cz(3);
    Matrix result(3);
    Matrix compare(3);  
    
    Cx(1) = 1.45;      Cx(2) = -4.32;       Cx(3) = -3.89;
    Cy(1) = 0.5;       Cy(2) = 2.45;        Cy(3) = 4.5;
    Cz(1) = -3.5;      Cz(2) = 3.8;         Cz(3) = -0.9;

    compare(1) = 2.226632074174593;   compare(2) = -1.451022829849176;    compare(3) = -1.337154847233099;

    result = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);

    _assert(m_equals(compare.data, result.data, 1, 3, 1e-14));

    return 0;
}

int EccAnom_01(){
    double M = 0.87; 
    double e = -3.63; 
    double compare = 0.188782546113942;
    double result = EccAnom(M, e);
    _assert(abs(compare - result)<1e-7);

    return 0;
}

int EccAnom_02(){
    double M = -1.36; 
    double e = 4.22; 
    double compare = 3.488484545359726;
    double result = EccAnom(M, e);
    _assert(abs(compare - result)<1e-7);

    return 0;
}

int Geodetic_01(){
    Matrix r(3);    r(1) = 1.1;   r(2) = 0.54;   r(3) = 1234;
    double lon;
    double lat;
    double h;

    double lon_comp = 0.456348468557109;
    double lat_comp = 1.570768524427010;
    double h_comp = -6.355517616575113e+06;
    
    Geodetic(r, lon, lat, h);
    _assert(abs(lon_comp - lon)<1e-7);
    _assert(abs(lat_comp - lat)<1e-7);
    _assert(abs(h_comp - h)<1e-7);

    return 0;
}

int Geodetic_02(){
    Matrix r(3);    r(1) = -0.65;   r(2) = 1.54;   r(3) = 5421;
    double lon;
    double lat;
    double h;

    double lon_comp = 1.970189346820603;
    double lat_comp = 1.570761691977959;
    double h_comp = -6.351330616563200e+06;
    
    Geodetic(r, lon, lat, h);
    _assert(abs(lon_comp - lon)<1e-7);
    _assert(abs(lat_comp - lat)<1e-7);
    _assert(abs(h_comp - h)<1e-7);

    return 0;
}

int IERS_01(){
    
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
    
    //Mjd_UTC debe aparecer en la cuarta columna de la matriz eopdata
    double Mjd_UTC = 37668;
    char interp = 'l';

    double x_pole_comp = -1.06654161707287e-07; 
    double y_pole_comp = 1.04865684037674e-06; 
    double UT1_UTC_comp = 0.0311435; 
    double LOD_comp = 0.001496; 
    double dpsi_comp = 3.09762005271316e-07;
    double deps_comp = 3.18716513961409e-08; 
    double dx_pole_comp = 0; 
    double dy_pole_comp = 0; 
    double TAI_UTC_comp = 2;

    IERS(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, eop, Mjd_UTC, interp);
    
    _assert(abs(x_pole_comp - x_pole)<1e-7);
    _assert(abs(y_pole_comp - y_pole)<1e-7);
    _assert(abs(UT1_UTC_comp - UT1_UTC)<1e-7);
    _assert(abs(LOD_comp - LOD)<1e-7);
    _assert(abs(dpsi_comp - dpsi)<1e-7);
    _assert(abs(deps_comp - deps)<1e-7);
    _assert(abs(dx_pole_comp - dx_pole)<1e-7);
    _assert(abs(dy_pole_comp - dy_pole)<1e-7);
    _assert(abs(TAI_UTC_comp - TAI_UTC)<1e-7);

    return 0;
}

int IERS_02(){
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
    
    //Mjd_UTC debe aparecer en la cuarta columna de la matriz eopdata
    double Mjd_UTC = 37732;

    double x_pole_comp = -5.32422384594492e-08; 
    double y_pole_comp = 1.25866358262296e-06; 
    double UT1_UTC_comp = 0.0184953; 
    double LOD_comp = 0.001999; 
    double dpsi_comp = 3.17354187517491e-07;
    double deps_comp = 1.94846618437923e-08; 
    double dx_pole_comp = 0; 
    double dy_pole_comp = 0; 
    double TAI_UTC_comp = 2;

    IERS(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, eop, Mjd_UTC);
    
    _assert(abs(x_pole_comp - x_pole)<1e-7);
    _assert(abs(y_pole_comp - y_pole)<1e-7);
    _assert(abs(UT1_UTC_comp - UT1_UTC)<1e-7);
    _assert(abs(LOD_comp - LOD)<1e-7);
    _assert(abs(dpsi_comp - dpsi)<1e-7);
    _assert(abs(deps_comp - deps)<1e-7);
    _assert(abs(dx_pole_comp - dx_pole)<1e-7);
    _assert(abs(dy_pole_comp - dy_pole)<1e-7);
    _assert(abs(TAI_UTC_comp - TAI_UTC)<1e-7);

    return 0;
}


int MeanObliquity_01(){
    double result;
    double Mjd_TT = 25463;
    double compare = 0.409254869415201;

    
    result = MeanObliquity(Mjd_TT);
    _assert(abs(compare - result)<1e-7);

    return 0;
}

int MeanObliquity_02(){
    double result;
    double Mjd_TT = -34.6;
    double compare = 0.409413285127502;

    
    result = MeanObliquity(Mjd_TT);
    _assert(abs(compare - result)<1e-7);

    return 0;
}

int Mjday_01(){
    int yr = 2024;
    int mon = 4;
    int day = 28;
    int hr = 19;
    int min = 39; 
    int sec = 10;
    double result;
    double compare = 6.042881886574067e+04;

    result = Mjday(yr, mon, day, hr, min, sec);
    _assert(abs(compare - result)<1e-7);

    return 0;
}

int Mjday_02(){
    int yr = 1998;
    int mon = 6;
    int day = 2;
    double result;
    double compare = 50966;

    result =  Mjday(yr, mon, day);
    _assert(abs(compare - result)<1e-7);

    return 0;
}

int Mjday_TDB_01(){
    double Mjd_TT = 45735489.25;
    double result;
    double compare = 4.573548924999999e+07;

    result =  Mjday_TDB(Mjd_TT);
    _assert(abs(compare - result)<1e-7);

    return 0;
}

int Mjday_TDB_02(){
    double Mjd_TT = 45735e6;
    double result;
    double compare = 4.573500000000000e+10;

    result =  Mjday_TDB(Mjd_TT);
    _assert(abs(compare - result)<1e-7);

    return 0;
}

int Position_01(){
    double lon = 243;
    double lat = 23;
    double h = 1345;
    Matrix result(3);
    Matrix compare(3);  
    compare(1) = 1.553565612896802;     compare(2) = 3.032594239753788;     compare(3) = -5.375212122319764;
    compare = compare*1.0e+06;

    result =  Position(lon, lat, h);
    _assert(m_equals(compare.data, result.data, 1, 3, 1e-8));

    return 0;
}

int Position_02(){
    double lon = 55;
    double lat = -65;
    double h = 520;
    Matrix result(3);
    Matrix compare(3);  
    compare(1) = -0.079566424227238;     compare(2) = 3.595056740585125;     compare(3) = -5.250780421296342;
    compare = compare*1.0e+06;

    result =  Position(lon, lat, h);
    _assert(m_equals(compare.data, result.data, 1, 3, 1e-8));

    return 0;
}

int NutAngles_01(){
    double Mjd_TT = 45.4;
    double dpsi;
    double deps;

    double dpsi_comp = 3.99424411030664e-05;
    double deps_comp = 3.61595902628267e-05;
    
    NutAngles(Mjd_TT, dpsi, deps);

    _assert(abs(dpsi_comp - dpsi)<1e-10);
    _assert(abs(deps_comp - deps)<1e-10);

    return 0;
}


int NutAngles_02(){
    double Mjd_TT = -9.51;
    double dpsi;
    double deps;

    double dpsi_comp = 2.66224415163593e-05;
    double deps_comp = 3.98255663249084e-05;
    
    NutAngles(Mjd_TT, dpsi, deps);

    _assert(abs(dpsi_comp - dpsi)<1e-10);
    _assert(abs(deps_comp - deps)<1e-10);

    return 0;
}

int EqnEquinox_01(){
    double Mjd_TT = 34.56;
    double result;
    double compare = 3.47470072006179e-05;
    
    result = EqnEquinox(Mjd_TT);

    _assert(abs(compare - result)<1e-7);

    return 0;
}


int EqnEquinox_02(){
    double Mjd_TT = -231.0;
    double result;
    double compare = 1.50737643141602e-05;
    
    result = EqnEquinox(Mjd_TT);

    _assert(abs(compare - result)<1e-7);

    return 0;
}


int NutMatrix_01(){
    double Mjd_TT = 4.854;
    Matrix result(3, 3);
    Matrix compare(3, 3);
    
    compare(1,1) = 0.999999999600563;           compare(1,2) = -2.59284527155697e-05;           compare(1,3) = -1.12512253584835e-05;        
    compare(2,1) = 2.59280240617824e-05;        compare(2,2) = 0.999999998938183;               compare(2,3) = -3.80968826992345e-05;         
    compare(3,1) = 1.12522131397585e-05;        compare(3,2) = 3.80965909620445e-05;            compare(3,3) = 0.999999999211019;
    
    result = NutMatrix(Mjd_TT);

    _assert(m_equals(compare.data, result.data, 1, 3, 1e-10));

    return 0;
}


int NutMatrix_02(){
    double Mjd_TT = -54.09;
    Matrix result(3, 3);
    Matrix compare(3, 3);
    
    compare(1,1) = 0.999999999562003;           compare(1,2) = -2.71511225037766e-05;           compare(1,3) = -1.17817946061955e-05;        
    compare(2,1) = 2.71506022766636e-05;        compare(2,2) = 0.999999998656675;               compare(2,3) = -4.41530810371127e-05;         
    compare(3,1) = 1.17829933960808e-05;        compare(3,2) = 4.41527611349524e-05;            compare(3,3) = 0.999999998955847;
    
    result = NutMatrix(Mjd_TT);

    _assert(m_equals(compare.data, result.data, 1, 3, 1e-10));

    return 0;
}


int PoleMatrix_01(){
    double xp = -54.09;
    double yp = 23.1;

    Matrix result(3, 3);
    Matrix compare(3, 3);
    
    compare(1,1) = -0.775730245535253;          compare(1,2) = -0.564921098476534;           compare(1,3) = -0.281259201907925;        
    compare(2,1) = 0;                           compare(2,2) = -0.445690000444333;           compare(2,3) = 0.895187367819682;         
    compare(3,1) = -0.631064644994328;          compare(3,2) = 0.694423916638818;            compare(3,3) = 0.345735213477289;
    
    result = PoleMatrix(xp, yp);

    _assert(m_equals(compare.data, result.data, 1, 3, 1e-10));

    return 0;
}


int PrecMatrix_01(){
    double Mjd_1 = 45.908;
    double Mjd_2 = 2.36;

    Matrix result(3, 3);
    Matrix compare(3, 3);
    
    compare(1,1) = 0.999999999577998;          compare(1,2) = 2.66386402804907e-05;        compare(1,3) = 1.15925392281205e-05;        
    compare(2,1) = -2.66386402804907e-05;      compare(2,2) = 0.999999999645191;           compare(2,3) = -1.54404709603792e-10;         
    compare(3,1) =  -1.15925392281205e-05;     compare(3,2) = -1.54404772896746e-10;       compare(3,3) = 0.999999999932807;
    
    result = PrecMatrix(Mjd_1, Mjd_2);

    _assert(m_equals(compare.data, result.data, 1, 3, 1e-10));

    return 0;
}


int TimeUpdate_01(){
    double Qdt = 45.908;

    Matrix P(3,3);
    P(1,1) = 1;         P(1,2) = 2;         P(1,3) = 3;         
    P(2,1) = 4;         P(2,2) = 5;         P(2,3) = 6;         
    P(3,1) = 7;         P(3,2) = 8;         P(3,3) = 9;   

    Matrix Phi(3,3);
    Phi(1,1) = 4;         Phi(1,2) = 0;         Phi(1,3) = 7.21;        
    Phi(2,1) = -2.3;      Phi(2,2) = -1.66;     Phi(2,3) = 0.32;         
    Phi(3,1) = 4.2;       Phi(3,2) = -8.61;     Phi(3,3) = -7.4;

    Matrix result(3, 3);
    Matrix compare(3, 3);
    
    compare(1,1) = 818.1649;          compare(1,2) = -163.797;          compare(1,3) = -859.8088;        
    compare(2,1) = -81.6878;          compare(2,2) = 74.0088;           compare(2,3) = 210.2474;         
    compare(3,1) = -1044.0606;        compare(3,2) = 356.5798;          compare(3,3) = 1291.2725;
    
    result = TimeUpdate(P, Phi, Qdt);

    _assert(m_equals(compare.data, result.data, 1, 3, 1e-10));

    return 0;
}

int TimeUpdate_02(){
    Matrix P(3,3);
    P(1,1) = -1;           P(1,2) = 2.1;           P(1,3) = 3;         
    P(2,1) = 0.4;          P(2,2) = -4.45;         P(2,3) = 6.54;         
    P(3,1) = 7.21;         P(3,2) = -2.8;          P(3,3) = -9.1;   

    Matrix Phi(3,3);
    Phi(1,1) = 4;         Phi(1,2) = 0;         Phi(1,3) = 7.21;        
    Phi(2,1) = -2.3;      Phi(2,2) = -1.66;     Phi(2,3) = 0.32;         
    Phi(3,1) = 4.2;       Phi(3,2) = -8.61;     Phi(3,3) = -7.4;

    Matrix result(3, 3);
    Matrix compare(3, 3);
    
    compare(1,1) = -194.59891;           compare(1,2) = -107.95087;           compare(1,3) = 699.7493;        
    compare(2,1) = -133.246364;          compare(2,2) = -18.440508;           compare(2,3) = 155.20639;         
    compare(3,1) = -73.615374;           compare(3,2) = 35.218722;            compare(3,3) = -1015.285285;
    
    result = TimeUpdate(P, Phi);

    _assert(m_equals(compare.data, result.data, 1, 3, 1e-10));

    return 0;
}

int angl_01(){
    double result;
    double compare = 2.63177277888791;

    Matrix vec1(3);     vec1(1)=-3.645;     vec1(2)=1.34;       vec1(3)=3.23;
    Matrix vec2(3);     vec2(1)=5.232;      vec2(2)=0.43;       vec2(3)=-1.89;

    result = angl(vec1, vec2);

    _assert(abs(compare - result)<1e-7);

    return 0;
}

int gmst_01(){
    double result;
    double compare = 3.47659871857155;
    double Mjd_UT1 = 21.34;

    result = gmst(Mjd_UT1);

    _assert(abs(compare - result)<1e-7);

    return 0;
}

int gmst_02(){
    double result;
    double compare = 0.425074383603393;
    double Mjd_UT1 = -0.087;

    result = gmst(Mjd_UT1);

    _assert(abs(compare - result)<1e-7);

    return 0;
}

int LTC_01(){
    Matrix result(3, 3);
    Matrix compare(3, 3);

    compare(1,1) = -0.606100368397834;          compare(1,2) = -0.795388171541424;           compare(1,3) = 0.0;        
    compare(2,1) = 0.784004541554557;           compare(2,2) = -0.597425833654158;           compare(2,3) = -0.168580105897647;         
    compare(3,1) = 0.134086622188189;           compare(3,2) = -0.10217646428911;            compare(3,3) = 0.985687956655421;

    double lon = 21.34;
    double lat = -4.543; 

    result = LTC(lon, lat);

    _assert(m_equals(compare.data, result.data, 1, 3, 1e-10));

    return 0;
}

int LTC_02(){
    Matrix result(3, 3);
    Matrix compare(3, 3);

    compare(1,1) = -0.931961590907104;           compare(1,2) = 0.362557020444922;             compare(1,3) = 0.0;        
    compare(2,1) = -0.328811740160542;           compare(2,2) = -0.84521853167509;             compare(2,3) = -0.421293808696023;         
    compare(3,1) = -0.152743028012723;           compare(3,2) = -0.392629648191659;            compare(3,3) = 0.906924212243999;

    double lon = -432.34;
    double lat = -54.543; 

    result = LTC(lon, lat);

    _assert(m_equals(compare.data, result.data, 1, 3, 1e-10));

    return 0;
}

int MeasUpdate_01(){
    Matrix K(3, 3);     
    Matrix x(3, 3);     
    x(1,1) = -3.21;     x(1,2) = -0.23;     x(1,3) = 2.13;    
    x(2,1) = 9.23;      x(2,2) = 1.2;       x(2,3) = 4.23;      
    x(3,1) = 1.23;      x(3,2) = 6.52;      x(3,3) = -3.01;

    Matrix P(3, 3);     
    P(1,1) = 8.32;      P(1,2) = 2.34;      P(1,3) = 2.36;
    P(2,1) = 2;         P(2,2) = 4.3;       P(2,3) = -2.14;
    P(3,1) = -5.3;      P(3,2) = -1.5;      P(3,3) = 0.28;

    double z = 2.34;

    double g = -4.52;

    double s = 1.48;

    Matrix G(3);     
    G(1) = -2.22;     G(2) = 3.12;      G(3) = -2.14;

    double n = 3;
    
    Matrix K_res(3);     
    K_res(1) = -0.243505079978358;      K_res(2) = 0.203505392241346;        K_res(3) = 0.097384016818965;      

    Matrix x_res(3, 3);     
    x_res(1,1) = -4.880444848651535;      x_res(1,2) = -1.900444848651534;     x_res(1,3) = 0.459555151348465;
    x_res(2,1) = 10.626046990775631;      x_res(2,2) = 2.596046990775631;      x_res(2,3) = 5.626046990775631;
    x_res(3,1) = 1.898054355378099;       x_res(3,2) = 7.188054355378099;      x_res(3,3) = -2.341945644621901;

    Matrix P_res(3, 3);    
    P_res(1,1) = 8.103670086947227;      P_res(1,2) = 5.123555270248607;      P_res(1,3) = -0.687514776945144;
    P_res(2,1) = 2.180794190467212;      P_res(2,2) = 1.973689160210730;      P_res(2,3) = 0.406910684978889;
    P_res(3,1) = -5.213484039458031;     P_res(3,2) = -2.613216173060952;     P_res(3,3) = 1.498780447292710;


    MeasUpdate(K, x, P, z, g, s, G, n);
    _assert(m_equals(K_res.data, K.data, 1, 3, 1e-10));
    _assert(m_equals(x_res.data, x.data, 1, 3, 1e-10));
    _assert(m_equals(P_res.data, P.data, 1, 3, 1e-10));

    return 0;
}

int elements_01(){
    
    Matrix y(6);    y(1) = 5.32;    y(2) = -0.54;   y(3) = -2.72;   y(4) = 3.95;    y(5) = -5.18;   y(6) = 0.58;
    double p;
    double a;
    double e;
    double i;
    double Omega;
    double omega;
    double M;

    double p_res = 2.62194587925332e-12;
    double a_res =  2.99968331661957;
    double e_res = 0.999999999999563;
    double i_res = 2.47584319520146;
    double Omega_res = 5.47748701838715;
    double omega_res = 2.31733513258329;
    double M_res = 3.14159136761281;


    elements(p, a, e, i, Omega, omega, M, y);

    _assert(abs(p_res - p)<1e-15);
    _assert(abs(a_res - a)<1e-6);
    _assert(abs(e_res - e)<1e-6);
    _assert(abs(i_res - i)<1e-6);
    _assert(abs(Omega_res - Omega)<1e-6);
    _assert(abs(omega_res - omega)<1e-6);
    _assert(abs(M_res - M)<1e-6);

    return 0;
}

int elements_02(){
    
    Matrix y(6);    y(1) = -3.62;    y(2) = 9.23;   y(3) = -1.45;   y(4) = 5.61;    y(5) = 4.98;   y(6) = -1.75;
    double p;
    double a;
    double e;
    double i;
    double Omega;
    double omega;
    double M;

    double p_res = 1.29510170235097e-11;
    double a_res =  5.00998502992149;
    double e_res = 0.999999999998707;
    double i_res = 2.90266198069576;
    double Omega_res = 5.73016838170846;
    double omega_res = 0.65790737182292;
    double M_res = 3.14159139173211;


    elements(p, a, e, i, Omega, omega, M, y);

    _assert(abs(p_res - p)<1e-15);
    _assert(abs(a_res - a)<1e-6);
    _assert(abs(e_res - e)<1e-6);
    _assert(abs(i_res - i)<1e-6);
    _assert(abs(Omega_res - Omega)<1e-6);
    _assert(abs(omega_res - omega)<1e-6);
    _assert(abs(M_res - M)<1e-6);

    return 0;
}


int AccelHarmonic_01(){
    Matrix result(3);
    Matrix r(3);
    r(1) = -3.21;     r(2) = -0.23;     r(3) = 2.13;

    Matrix E(3, 3);     
    E(1,1) = 8.32;      E(1,2) = 2.34;      E(1,3) = 2.36;
    E(2,1) = 2;         E(2,2) = 4.3;       E(2,3) = -2.14;
    E(3,1) = -5.3;      E(3,2) = -1.5;      E(3,3) = 0.28;

    Matrix compare(3);
    compare(1) = -6.49172129990248e+23;      compare(2) = -2.69320945943847e+23;      compare(3) = -8.20812319455523e+22;

    double n_max = 3;
    double m_max = 3;

    result = AccelHarmonic(r, E, n_max, m_max);

    _assert(m_equals(compare.data, result.data, 1, 3, 1e10));

    return 0;
}

int AccelHarmonic_02(){
    Matrix result(3);
    Matrix r(3);
    r(1) = 34.64;     r(2) = -76.31;     r(3) =65.25;
    

    Matrix E(3, 3);     
    E(1,1) = 145.24;       E(1,2) = 16.54;       E(1,3) = 37.76;
    E(2,1) = 74.67;        E(2,2) = 28.34;       E(2,3) = 28.25;
    E(3,1) = -72.56;       E(3,2) = 86.47;       E(3,3) = -65.13;

    Matrix compare(3);
    compare(1) = 51271633627.927;      compare(2) = 177500101964.719;      compare(3) = -58735951566.3968;

    double n_max = 3;
    double m_max = 3;

    result = AccelHarmonic(r, E, n_max, m_max);

    _assert(m_equals(compare.data, result.data, 1, 3, 1e3));

    return 0;
}

int AccelHarmonic_03(){
    Matrix result(3);
    Matrix r(3);
    r(1) = 4564.956;     r(2) = -4786.587;     r(3) = 6869.5684;
    
    Matrix E(3, 3);     
    E(1,1) = 13354.465;       E(1,2) = 8762.2414;       E(1,3) = -46548.2;
    E(2,1) = 65658.524;       E(2,2) = -4542.54;        E (2,3) = 984.54;
    E(3,1) = 5968.45;         E(3,2) = -4556.756;       E(3,3) = 4873.8564;

    Matrix compare(3);
    compare(1) = -77.4203028224361;      compare(2) = 19.3352573942419;      compare(3) = -63.2347204300577;

    double n_max = 20;
    double m_max = 20;

    result = AccelHarmonic(r, E, n_max, m_max);

    _assert(m_equals(compare.data, result.data, 1, 3, 1e-9));

    return 0;
}

int gast_01(){
    double Mjd_UT1 = 234.64;
    double compare = 2.747772712783872;
    double result = gast(Mjd_UT1);

    _assert(abs(compare - result)<1e-6);

    return 0;
}

int GHAMatrix_01(){
    Matrix result(3, 3);         
    double Mjd_UT1 = 156.265;

    Matrix compare(3, 3);     
    compare(1,1) = 0.576218856780285;       compare(1,2) = -0.817295435623387;      compare(1,3) = 0.0;
    compare(2,1) = 0.817295435623387;       compare(2,2) = 0.576218856780285;       compare(2,3) = 0.0;
    compare(3,1) = 0.0;                     compare(3,2) = 0.0;                     compare(3,3) = 1.0;   
    result = GHAMatrix(Mjd_UT1);

    _assert(m_equals(compare.data, result.data, 3, 3, 1e-6));

    return 0;
}

int GHAMatrix_02(){
    Matrix result(3, 3);
    double Mjd_UT1 = 49746.1101542334;

    Matrix compare(3, 3);
    compare(1,1) = -0.976451404712871;       compare(1,2) = 0.215737465995733;      compare(1,3) = 0.0;
    compare(2,1) = -0.215737465995733;       compare(2,2) = -0.976451404712871;       compare(2,3) = 0.0;
    compare(3,1) = 0.0;                     compare(3,2) = 0.0;                     compare(3,3) = 1.0;
    result = GHAMatrix(Mjd_UT1);

    _assert(m_equals(compare.data, result.data, 3, 3, 1e-6));

    return 0;
}

int G_AccelHarmonic_01(){
    Matrix result(3, 3);     
    Matrix r(3);
    r(1) = 1.43;     r(2) = -2.76;     r(3) = 0.54;
    
    Matrix E(3, 3);     
    E(1,1) = 23.43;      E(1,2) = 72.41;      E(1,3) = 47.13;
    E(2,1) = -13.54;     E(2,2) = 13.67;      E(2,3) = -25.76;
    E(3,1) = 24.67;      E(3,2) = 82.25;      E(3,3) = 26.14;

    Matrix compare(3, 3);     
    compare(1,1) = 1.6145357703404e+19;         compare(1,2) = -4.4410883933236e+19;      compare(1,3) = 3.58877206425741e+19;
    compare(2,1) = -4.83719961813644e+19;       compare(2,2) = -5.94041543085921e+20;      compare(2,3) = -3.86395883692299e+19;
    compare(3,1) = 3.16652795172409e+19;        compare(3,2) = -3.85802453564528e+19;      compare(3,3) = 7.15842564866154e+19;      

    double n_max = 3;
    double m_max = 3;

    result = G_AccelHarmonic(r, E, n_max, m_max);

    _assert(m_equals(compare.data, result.data, 3, 3, 1e8));

    return 0;
}


int G_AccelHarmonic_02(){
    Matrix result(3, 3);     
    Matrix r(3);
    r(1) = 0.25;     r(2) = 0.63;     r(3) = -1.23;
    
    Matrix E(3, 3);     
    E(1,1) = 2.13;      E(1,2) = 1.36;      E(1,3) =-0.55;
    E(2,1) = 0.28;      E(2,2) = -2.92;     E(2,3) = -1.28;
    E(3,1) = -1.64;     E(3,2) = 2.19;      E(3,3) = 0.34;

    Matrix compare(3, 3);     
    compare(1,1) = 6.83981855530016e+28;        compare(1,2) = -2.40263603532317e+28;      compare(1,3) = -1.5298291616873e+28;
    compare(2,1) = 5.54265849226182e+29;        compare(2,2) = -5.77105436729858e+27;      compare(2,3) = 8.94866801436029e+28;
    compare(3,1) = 6.94094619842872e+28;        compare(3,2) = -3.5923441629265e+28;       compare(3,3) = -3.63506461770632e+27;      

    double n_max = 3;
    double m_max = 3;

    result = G_AccelHarmonic(r, E, n_max, m_max);

    _assert(m_equals(compare.data, result.data, 3, 3, 1e18));

    return 0;
}

int JPL_Eph_DE430_01(){

    Matrix aux = (*Global::PC);
    // printf("%5.20lf",aux(1,1));
    // printf("%5.20lf",aux(1,2));
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
    double Mjd_TDB = 37756;

    Matrix r_Mercury_comp(3);   r_Mercury_comp(1) = 193857295552.949;        r_Mercury_comp(2) = -1586683282.33992;       r_Mercury_comp(3) = -8742151553.24134;
    Matrix r_Venus_comp(3);     r_Venus_comp(1) = 215165189050.982;          r_Venus_comp(2) = 106628724393.04;           r_Venus_comp(3) = 43114261974.4799;
    Matrix r_Earth_comp(3);     r_Earth_comp(1) = -146652871565.183;         r_Earth_comp(2) = -28458528909.0728;         r_Earth_comp(3) = -12335597391.8487;
    Matrix r_Mars_comp(3);      r_Mars_comp(1) = 322423116458.008;           r_Mars_comp(2) = -67071731672.6544;          r_Mars_comp(3) = -36275010440.0661;
    Matrix r_Jupiter_comp(3);   r_Jupiter_comp(1) = 760167232265.62;         r_Jupiter_comp(2) = -365258705923.513;       r_Jupiter_comp(3) = -171399631149.296;
    Matrix r_Saturn_comp(3);    r_Saturn_comp(1) = 999899658456.895;         r_Saturn_comp(2) = -1089146069379.98;        r_Saturn_comp(3) = -485813405409.846;
    Matrix r_Uranus_comp(3);    r_Uranus_comp(1) = -2218955556426.83;        r_Uranus_comp(2) = 1288952894493.1;          r_Uranus_comp(3) = 597924262625.312;
    Matrix r_Neptune_comp(3);   r_Neptune_comp(1) = -3198932852767.42;       r_Neptune_comp(2) = -2835721760900.98;       r_Neptune_comp(3) = -1076733105159.76;
    Matrix r_Pluto_comp(3);     r_Pluto_comp(1) = -4399671140908.45;         r_Pluto_comp(2) = 1132898893472.76;          r_Pluto_comp(3) = 1726621192596.08;
    Matrix r_Moon_comp(3);      r_Moon_comp(1) = 320845955.347685;           r_Moon_comp(2) = -149448843.685816;          r_Moon_comp(3) = -74166152.0408648;
    Matrix r_Sun_comp(3);       r_Sun_comp(1) = 146099041202.014;            r_Sun_comp(2) = 29246385387.3131;            r_Sun_comp(3) = 12683942495.7121;
    
    JPL_Eph_DE430(r_Mercury, r_Venus, r_Earth, r_Mars,
                   r_Jupiter, r_Saturn, r_Uranus, r_Neptune,
                   r_Pluto, r_Moon, r_Sun, Mjd_TDB);


    _assert(m_equals(r_Mercury_comp.data, r_Mercury.data, 1, 3, 1e3));
    _assert(m_equals(r_Venus_comp.data, r_Venus.data, 1, 3, 1e3));
    _assert(m_equals(r_Earth_comp.data, r_Earth.data, 1, 3, 1e3));
    _assert(m_equals(r_Mars_comp.data, r_Mars.data, 1, 3, 1e14));
    _assert(m_equals(r_Jupiter_comp.data, r_Jupiter.data, 1, 3, 1e3));
    _assert(m_equals(r_Saturn_comp.data, r_Saturn.data, 1, 3, 1e3));
    _assert(m_equals(r_Uranus_comp.data, r_Uranus.data, 1, 3, 1e3));
    _assert(m_equals(r_Neptune_comp.data, r_Neptune.data, 1, 3, 1e3));
    _assert(m_equals(r_Pluto_comp.data, r_Pluto.data, 1, 3, 1e3));
    _assert(m_equals(r_Moon_comp.data, r_Moon.data, 1, 3, 1e3));
    _assert(m_equals(r_Sun_comp.data, r_Sun.data, 1, 3, 1e3));

    return 0;
}

int Accel_01(){
    Global::auxparam.Mjd_UTC=49746.1163541665;

    double Mdj_TT = 37691.0;
    Matrix Y(6);
    Y(1) = 235654.1566;        Y(2) = 566544.4185;       Y(3) = 1546545.541;
    Y(4) = 954645.456;         Y(5) = 1385415.796;       Y(6) = -453368.454;

    Matrix result = Accel(Mdj_TT, Y);

    Matrix compare (6);
    compare(1) = 954645.456;               compare(2) = 1385415.796;            compare(3) = -453368.454;
    compare(4) = -12198826.5349416;        compare(5) =  -180686329.73599;      compare(6) = -165289008.651588;

    _assert(m_equals(compare.data, result.data, 1, 6, 1e2));

    return 0;
}

int gibbs_01(){
    Matrix v2(3);
    double theta, theta1, cop;
    string error = "";
    Matrix r1(3);       r1(1) = 5.157;       r1(2) = -1.276;    r1(3) = -3.573;
    Matrix r2(3);       r2(1) = -1.296;      r2(2) = 3.487;     r2(3) = 0.16;
    Matrix r3(3);       r3(1) = -1.394;      r3(2) = 0.225;     r3(3) = 4.133;

    Matrix v2_comp(3);       v2_comp(1) = -4670901.09790465;      v2_comp(2) = -2883373.07788591;     v2_comp(3) = 10138076.6027918;
    double theta_comp = 2.08401928517999;
    double theta1_comp = 1.36944044290861;
    double cop_comp = 0.526462288353585;
    string error_comp = "not coplanar";

    gibbs(r1, r2, r3, v2, theta, theta1, cop, error);

    _assert(m_equals(v2_comp.data, v2.data, 1, 3, 1e-3));
    _assert(abs(theta_comp - theta) < 1e-9);
    _assert(abs(theta1_comp - theta1) < 1e-9);
    _assert(abs(cop_comp - cop) < 1e-9);
    _assert(error.compare(error_comp) == 0);
    return 0;
}

int gibbs_02(){
    Matrix v2(3);
    double theta, theta1, cop;
    string error = "";
    Matrix r1(3);       r1(1) = 6.658;      r1(2) = 5.349;    r1(3) = 6.364;
    Matrix r2(3);       r2(1) = 3.248;      r2(2) = 8.33;     r2(3) = 5.546;
    Matrix r3(3);       r3(1) = 3.395;      r3(2) = 4.648;    r3(3) = 3.921;

    Matrix v2_comp(3);       v2_comp(1) = 0;      v2_comp(2) = 0;     v2_comp(3) = 0;
    double theta_comp = 0;
    double theta1_comp = 0;
    double cop_comp = -0.0319785534944269;
    string error_comp = "impossible";

    gibbs(r1, r2, r3,v2, theta, theta1, cop, error);

    _assert(m_equals(v2_comp.data, v2.data, 1, 3, 1e-3));
    _assert(abs(theta_comp - theta) < 1e-9);
    _assert(abs(theta1_comp - theta1) < 1e-9);
    _assert(abs(cop_comp - cop) < 1e-9);
    _assert(error.compare(error_comp) == 0);
    return 0;
}

int hgibbs_01(){
    Matrix v2(3);
    double theta, theta1, cop;
    string error = "";

    Matrix r1(3);       r1(1) = 6.658;      r1(2) = 5.349;    r1(3) = 6.364;
    Matrix r2(3);       r2(1) = 3.248;      r2(2) = 8.33;     r2(3) = 5.546;
    Matrix r3(3);       r3(1) = 3.395;      r3(2) = 4.648;    r3(3) = 3.921;
    double Mjd1 = 37675.0;
    double Mjd2 = 37677.0;
    double Mjd3 = 37693.0;

    Matrix v2_comp(3);       v2_comp(1) = -8.32966565160303e+16;      v2_comp(2) = 1.63061646963647e+17;     v2_comp(3) = 1.60911896731326e+16;
    double theta_comp = 0.438106733276759;
    double theta1_comp = 0.221143618458218;
    double cop_comp = -0.0319785534944269;
    string error_comp = "angl > 1ø";

    hgibbs(r1, r2, r3, Mjd1, Mjd2, Mjd3, v2, theta, theta1, cop, error);

    _assert(m_equals(v2_comp.data, v2.data, 1, 3, 1e8));
    _assert(abs(theta_comp - theta) < 1e-9);
    _assert(abs(theta1_comp - theta1) < 1e-9);
    _assert(abs(cop_comp - cop) < 1e-9);
    _assert(error.compare(error_comp) == 0);
    return 0;
}


int VarEqn_01(){
    double x = 42.015;
    Matrix yPhi(42);
    yPhi(1) = 7101576.98990384;       yPhi(2) = 1295199.87127754;       yPhi(3) = 12739.2823333892;       yPhi(4) = 576.004651193009;
    yPhi(5) = -3084.62203617271;      yPhi(6) = -6736.02594582756;      yPhi(7) = 1.00002525535511;       yPhi(8) = 7.08259815373561e-06;
    yPhi(9) = 1.91608861002907e-07;   yPhi(10) = 1.01043851887223e-05;  yPhi(11) = 2.82768336557965e-06;  yPhi(12) = 6.44131451075285e-08;
    yPhi(13) = 7.08259834024473e-06;  yPhi(14) = 0.999988040046622;     yPhi(15) = 3.53015288644891e-08;  yPhi(16) = 2.82768357826951e-06;
    yPhi(17) = -4.78603729288896e-06; yPhi(18) = 1.18527461137171e-08;  yPhi(19) = 1.9160935046062e-07;   yPhi(20) = 3.53016114843062e-08;
    yPhi(21) = 0.999986704774626;     yPhi(22) = 6.44136325079115e-08;  yPhi(23) = 1.18528331537947e-08;  yPhi(24) = -5.31820682446032e-06;
    yPhi(25) = 5.00001498082565;      yPhi(26) = 1.1781862826826e-05;   yPhi(27) = 2.68389762645616e-07;  yPhi(28) = 1.00002526606744;
    yPhi(29) = 7.05571100144785e-06;  yPhi(30) = 1.30455137405173e-07;  yPhi(31) = 1.17818628919961e-05;  yPhi(32) = 4.99995293819715;
    yPhi(33) = 4.93630678596596e-08;  yPhi(34) = 7.05571117883108e-06;  yPhi(35) = 0.999988029832331;     yPhi(36) = 2.39618837211068e-08;
    yPhi(37) = 2.68390168073246e-07;  yPhi(38) = 4.93631303180711e-08;  yPhi(39) = 4.99995072081276;      yPhi(40) = 1.30455621823661e-07;
    yPhi(41) = 2.39619698989173e-08;  yPhi(42) = 0.999986704276552;
    
    Matrix compare(42);
    compare(1) = 576.004651193009;          compare(2) = -3084.62203617271;         compare(3) = -6736.02594582756;         compare(4) = -7.53466223591456;
    compare(5) = -1.37422019436676;         compare(6) = -0.0135523187845607;       compare(7) = 1.01043851887223e-05;      compare(8) = 2.82768336557965e-06;
    compare(9) = 6.44131451075285e-08;      compare(10) = 2.02219654159752e-06;     compare(11) = 5.62315203773653e-07;     compare(12) = 5.54306320677402e-09;
    compare(13) = 2.82768357826951e-06;     compare(14) = -4.78603729288896e-06;    compare(15) = 1.18527461137171e-08;     compare(16) = 5.62315388477388e-07;
    compare(17) = -9.58426101464243e-07;    compare(18) = 1.01508482681555e-09;     compare(19) = 6.44136325079115e-08;     compare(20) = 1.18528331537947e-08;
    compare(21) = -5.31820682446032e-06;    compare(22) = 5.54345424372027e-09;     compare(23) = 1.01515597613669e-09;     compare(24) = -1.06368579634054e-06;
    compare(25) = 1.00002526606744;         compare(26) = 7.05571100144785e-06;     compare(27) = 1.30455137405173e-07;     compare(28) = 1.01107443600397e-05;
    compare(29) = 2.81153608464334e-06;     compare(30) = 2.77154087047314e-08;     compare(31) = 7.05571117883108e-06;     compare(32) = 0.999988029832331;
    compare(33) = 2.39618837211068e-08;     compare(34) = 2.81153631885262e-06;     compare(35) = -4.79215600401207e-06;    compare(36) = 5.07544131782405e-09;
    compare(37) = 1.30455621823661e-07;     compare(38) = 2.39619698989173e-08;     compare(39) = 0.999986704276552;        compare(40) = 2.77159004655376e-08;
    compare(41) = 5.07553139893108e-09;     compare(42) = -5.31844727803824e-06;
    
    Matrix result = VarEqn(x, yPhi);

    _assert(m_equals(compare.data, result.data, 1, 42, 1e-3));

    return 0;
}

int DEInteg_01(){
    int neqn = 6;
    Matrix y (6);
    y(1) = 6221397.62857869;    y(2) = 2867713.77965738;    y(3) = 3006155.98509949;
    y(4) = 4645.04725161806;    y(5) = -2752.21591588204;   y(6) = -7507.99940987031;
    double t = 0;
    double tout = -134.999991953373;
    double relerr = 1e-13;
    double abserr = 1e-6;

    int iflag = 1;
    double *work = new double[100+21*neqn];
    int iwork[5];

    Matrix y_comp(6);
    y_comp(1) = 5542555.89427452;    y_comp(2) = 3213514.83814162;    y_comp(3) = 3990892.92789074;
    y_comp(4) = 5394.06894044389;    y_comp(5) = -2365.2129057402;   y_comp(6) = -7061.8448137347;

    DEInteg (Accel, neqn, y, t, tout, relerr, abserr, iflag, work, iwork);

    _assert(m_equals(y.data, y_comp.data, 1, 6, 1e-1));

    return 0;
}

int DEInteg_02(){
    int neqn = 42;
    Matrix y (42);
    y(1) = 7101576.98990384;     y(2) = 1295199.87127754;        y(3) = 12739.2823333892;        y(4) = 576.004651193009;
    y(5) = -3084.62203617271;       y(6) = -6736.02594582756;       y(7) = 1.00002525535511;        y(8) = 7.08259815373561e-06;
    y(9) = 1.91608861002907e-07;    y(10) = 1.01043851887223e-05;   y(11) = 2.82768336557965e-06;   y(12) = 6.44131451075285e-08;
    y(13) = 7.08259834024473e-06;   y(14) = 0.999988040046622;      y(15) = 3.53015288644891e-08;   y(16) = 2.82768357826951e-06;
    y(17) = -4.78603729288896e-06;  y(18) = 1.18527461137171e-08;   y(19) = 1.9160935046062e-07;    y(20) = 3.53016114843062e-08;
    y(21) = 0.999986704774626;      y(22) = 6.44136325079115e-08;   y(23) = 1.18528331537947e-08;   y(24) = -5.31820682446032e-06;
    y(25) = 5.00001498082565;       y(26) = 1.1781862826826e-05;    y(27) = 2.68389762645616e-07;   y(28) = 1.00002526606744;
    y(29) = 7.05571100144785e-06;   y(30) = 1.30455137405173e-07;   y(31) = 1.17818628919961e-05;   y(32) = 4.99995293819715;
    y(33) = 4.93630678596596e-08;   y(34) = 7.05571117883108e-06;   y(35) = 0.999988029832331;      y(36) = 2.39618837211068e-08;
    y(37) = 2.68390168073246e-07;   y(38) = 4.93631303180711e-08;   y(39) = 4.99995072081276;       y(40) = 1.30455621823661e-07;
    y(41) = 2.39619698989173e-08;   y(42) = 0.999986704276552;

    double t = 0;
    double tout = 4.99997287988663;
    double relerr = 1e-13;
    double abserr = 1e-6;

    int iflag = 1;
    double *work = new double[100+21*neqn];
    int iwork[5];

    Matrix y_comp(42);
    y_comp(1) = 7104362.80284458;     y_comp(2) = 1279759.73542038;         y_comp(3) = -20940.6848441728;       y_comp(4) = 538.324123651396;
    y_comp(5) = -3091.45215137201;       y_comp(6) = -6736.00414345875;        y_comp(7) = 1.00010106496325;        y_comp(8) = 2.82228797112616e-05;
    y_comp(9) = 5.21764823246293e-07;    y_comp(10) = 2.02217660695483e-05;    y_comp(11) = 5.62308617974723e-06;   y_comp(12) = 5.54079508472868e-08;
    y_comp(13) = 2.82228869257608e-05;   y_comp(14) = 0.999952119814454;       y_comp(15) = 9.61122253578372e-08;   y_comp(16) = 5.62308983126859e-06;
    y_comp(17) = -9.58411796113994e-06;  y_comp(18) = 1.02566514140829e-08;    y_comp(19) = 5.21780417252776e-07;   y_comp(20) = 9.61150455025847e-08;
    y_comp(21) = 0.999946818043692;      y_comp(22) = 5.54157547633479e-08;    y_comp(23) = 1.02580674631061e-08;   y_comp(24) = -1.06365195667e-05;
    y_comp(25) = 10.0002827788113;       y_comp(26) = 9.37163424370079e-05;    y_comp(27) = 9.23630171094507e-07;   y_comp(28) = 1.00010114792979;
    y_comp(29) = 2.80071897080942e-05;   y_comp(30) = 3.23301815979299e-08;    y_comp(31) = 9.37163542143084e-05;   y_comp(32) = 9.9997860235401;
    y_comp(33) = 1.70239737106332e-07;   y_comp(34) = 2.80071969374962e-05;    y_comp(35) = 0.999952038496836;      y_comp(36) = 6.45715077098127e-09;
    y_comp(37) = 9.23656078243316e-07;   y_comp(38) = 1.70244437291363e-07;    y_comp(39) = 9.99976848260984;       y_comp(40) = 3.2345797180921e-08;
    y_comp(41) = 6.45998875635734e-09;   y_comp(42) = 0.999946816394809;

    DEInteg (VarEqn, neqn, y, t, tout, relerr, abserr, iflag, work, iwork);

    _assert(m_equals(y.data, y_comp.data, 1, 42, 1e-1));

    return 0;
}

int anglesg_01(){
    double az1 = (*Global::obs)(1,2);
    double az2 = (*Global::obs)(9,2);
    double az3 = (*Global::obs)(18,2);
    double el1 = (*Global::obs)(1,3);
    double el2 = (*Global::obs)(9,3);
    double el3 = (*Global::obs)(18,3);
    double Mjd1 = (*Global::obs)(1,1);
    double Mjd2 = (*Global::obs)(9,1);
    double Mjd3 = (*Global::obs)(18,1);

    Matrix Rs(3);   Rs(1) = -5512567.84003607;  Rs(2) = -2196994.44666933;  Rs(3) = 2330804.96614689;
    Matrix r2(3);
    Matrix v2(3);

    Matrix r2_comp(3);  r2_comp(1) = 6221397.62857869;  r2_comp(2) = 2867713.77965738;  r2_comp(3) = 3006155.98509949;
    Matrix v2_comp(3);  v2_comp(1) = 4645.04725161806;  v2_comp(2) = -2752.21591588204; v2_comp(3) = -7507.99940987031;

    anglesg (az1, az2, az3, el1, el2, el3, Mjd1, Mjd2, Mjd3, Rs, Rs, Rs, r2, v2);

    _assert(m_equals(r2.data, r2_comp.data, 1, 3, 1e1));
    _assert(m_equals(v2.data, v2_comp.data, 1, 3, 1e1));

    return 0;
}

int all_tests(){
    cout<< "Entro en r_x_01\n";
    _verify(r_x_01);
    cout<< "Entro en r_x_02\n";
    _verify(r_x_02);
    cout<< "Entro en r_y_01\n";
    _verify(r_y_01);
    cout<< "Entro en r_y_02\n";
    _verify(r_y_02);
    cout<< "Entro en r_z_01\n";
    _verify(r_z_01);
    cout<< "Entro en r_z_02\n";
    _verify(r_z_02);
    cout<< "Entro en legendre_01\n";
    _verify(legendre_01);
    cout<< "Entro en legendre_02\n";
    _verify(legendre_02);
    cout<< "Entro en legendre_03\n";
    _verify(legendre_03);
    cout<< "Entro en timediff_01\n";
    _verify(timediff_01);
    cout<< "Entro en timediff_02\n";
    _verify(timediff_02);
    cout<< "Entro en unit_01\n";
    _verify(unit_01);
    cout<< "Entro en unit_02\n";
    _verify(unit_02);
    cout<< "Entro en accelPointMass_01\n";
    _verify(accelPointMass_01);
    cout<< "Entro en accelPointMass_02\n";
    _verify(accelPointMass_02);
    cout<< "Entro en AzElPa_01\n";
    _verify(AzElPa_01);
    cout<< "Entro en AzElPa_02\n";
    _verify(AzElPa_02);
    cout<< "Entro en Cheb3D_01\n";
    _verify(Cheb3D_01);
    cout<< "Entro en Cheb3D_02\n";
    _verify(Cheb3D_02);
    cout<< "Entro en EccAnom_01\n";
    _verify(EccAnom_01);
    cout<< "Entro en EccAnom_02\n";
    _verify(EccAnom_02);
    cout<< "Entro en Geodetic_01\n";
    _verify(Geodetic_01);
    cout<< "Entro en Geodetic_02\n";
    _verify(Geodetic_02);
    cout<< "Entro en IERS_01\n";
    _verify(IERS_01);
    cout<< "Entro en IERS_02\n";
    _verify(IERS_02);
    cout<< "Entro en MeanObliquity_01\n";
    _verify(MeanObliquity_01);
    cout<< "Entro en MeanObliquity_02\n";
    _verify(MeanObliquity_02);
    cout<< "Entro en Mjday_01\n";
    _verify(Mjday_01);
    cout<< "Entro en Mjday_02\n";
    _verify(Mjday_02);
    cout<< "Entro en Mjday_TDB_01\n";
    _verify(Mjday_TDB_01);
    cout<< "Entro en Mjday_TDB_02\n";
    _verify(Mjday_TDB_02);
    cout<< "Entro en Position_01\n";
    _verify(Position_01);
    cout<< "Entro en Position_02\n";
    _verify(Position_02);
    cout<< "Entro en NutAngles_01\n";
    _verify(NutAngles_01);
    cout<< "Entro en NutAngles_02\n";
    _verify(NutAngles_02);
    cout<< "Entro en EqnEquinox_01\n";
    _verify(EqnEquinox_01);
    cout<< "Entro en EqnEquinox_02\n";
    _verify(EqnEquinox_02);
    cout<< "Entro en NutMatrix_01\n";
    _verify(NutMatrix_01);
    cout<< "Entro en NutMatrix_02\n";
    _verify(NutMatrix_02);
    cout<< "Entro en PoleMatrix_01\n";
    _verify(PoleMatrix_01);
    cout<< "Entro en PrecMatrix_01\n";
    _verify(PrecMatrix_01);
    cout<< "Entro en TimeUpdate_01\n";
    _verify(TimeUpdate_01);
    cout<< "Entro en TimeUpdate_02\n";
    _verify(TimeUpdate_02);
    cout<< "Entro en angl_01\n";
    _verify(angl_01);
    cout<< "Entro en gmst_01\n";
    _verify(gmst_01);
    cout<< "Entro en gmst_02\n";
    _verify(gmst_02);
    cout<< "Entro en LTC_01\n";
    _verify(LTC_01);
    cout<< "Entro en LTC_02\n";
    _verify(LTC_02);
    cout<< "Entro en MeasUpdate_01\n";
    _verify(MeasUpdate_01);
    cout<< "Entro en elements_01\n";
    _verify(elements_01);
    cout<< "Entro en elements_02\n";
    _verify(elements_02);
    cout<< "Entro en AccelHarmonic_01\n";
    _verify(AccelHarmonic_01);
    cout<< "Entro en AccelHarmonic_02\n";
    _verify(AccelHarmonic_02);
    cout<< "Entro en AccelHarmonic_03\n";
    _verify(AccelHarmonic_03);
    cout<< "Entro en gast_01\n";
    _verify(gast_01);
    cout<< "Entro en GHAMatrix_01\n";
    _verify(GHAMatrix_01);
    cout<< "Entro en GHAMatrix_02\n";
    _verify(GHAMatrix_02);
    cout<< "Entro en G_AccelHarmonic_01\n";
    _verify(G_AccelHarmonic_01);
    cout<< "Entro en G_AccelHarmonic_02\n";
    _verify(G_AccelHarmonic_02);
    cout<< "Entro en JPL_Eph_DE430_01\n";
    _verify(JPL_Eph_DE430_01);
    cout<< "Entro en Accel_01\n";
    _verify(Accel_01);
    cout<< "Entro en gibbs_01\n";
    _verify(gibbs_01);
    cout<< "Entro en gibbs_02\n";
    _verify(gibbs_02);
    cout<< "Entro en hgibbs_01\n";
    _verify(hgibbs_01);
    cout<< "Entro en VarEqn_01\n";
    _verify(VarEqn_01);
    cout<< "Entro en DEInteg_01\n";
    _verify(DEInteg_01);
    cout<< "Entro en DEInteg_02\n";
    _verify(DEInteg_02);
    cout<< "Entro en anglesg_01\n";
    _verify(anglesg_01);
    cout<<"Todos los test son correctos\n";

    return 0;
}
/*
int main() {
    Global::eop19620101(21413);
    Global::DE430Coeff(2285, 1020);
    Global::GEOS3(46);
    Global::GGM03S();
    double Mjd_UTC = (*Global::obs)(9, 1);
    Global::auxparam.Mjd_UTC = Mjd_UTC;
    Global::auxparam.n = 20;
    Global::auxparam.m = 20;
    Global::auxparam.sun = 1;
    Global::auxparam.moon = 1;
    Global::auxparam.planets = 1;
    all_tests();
    return 0;
}
*/
//Accel
//AuxParam es un conjunto de datos que se genera en el principal. Es un struct.
/*

*/

/*DEInteg
traga funciones como Global::auxparam
mirar archivos subidos en el aula
https://people.sc.fsu.edu/~jburkardt/f77_src/ode/ode.html

 int iflag;
 iflag = 1;
 double *work;
 work = new double[100+21*neqn];
 int iwork[5];
*/

/*anglesg
roots devuelve las raices de un polinomio
solucion: robar codigo de crbond.com/downloads/rpoly.cpp
void roots(Matrix coef, Matrix & real, Matrix & im){
    sacar vectores de las matrices
    llamar rpoly
    coger resutlados y adaptar a matrices
}
*/