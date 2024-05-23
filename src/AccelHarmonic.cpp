#include "..\include\AccelHarmonic.h"

Matrix AccelHarmonic(Matrix r, Matrix E, int n_max, int m_max){

    double r_ref = 6378.1363e3;   // Earth's radius [m]; GGM03S
    double gm    = 398600.4415e9; // [m^3/s^2]; GGM03S
    
    // Body-fixed position 
    Matrix r_bf(3);
    Matrix Maux = E * transpose(r);
    r_bf(1) = Maux(1,1);     r_bf(2) = Maux(2,1);     r_bf(3) = Maux(3,1);

    // Auxiliary quantities
    double d = norm(r_bf);  // distance
    double latgc = asin(r_bf(3)/d); 
    double lon = atan2(r_bf(2),r_bf(1));

    int aux;
    if (n_max>m_max)
        aux=n_max;
    else
        aux=m_max;
    Matrix pnm(n_max+1,aux+1);
    Matrix dpnm(n_max+1,aux+1);
    
    Legendre(n_max, m_max, latgc, pnm, dpnm);

    double dUdr = 0;
    double dUdlatgc = 0;
    double dUdlon = 0;
    double q3 = 0; double q2 = q3; double q1 = q2;
    double b1; double b2; double b3;
    
    for (int n=0; n<=n_max; n++){
        b1 = -gm/pow(d,2)*pow(r_ref/d, n)*(n + 1);
        b2 = (gm/d)*pow(r_ref/d, n);
        b3 = (gm/d)*pow(r_ref/d, n);
        for (int m=0; m<=m_max; m++){
            q1 = q1 + pnm(n+1,m+1)*((*Global::Cnm)(n+1,m+1)*cos(m*lon)+(*Global::Snm)(n+1,m+1)*sin(m*lon));
            q2 = q2 + dpnm(n+1,m+1)*((*Global::Cnm)(n+1,m+1)*cos(m*lon)+(*Global::Snm)(n+1,m+1)*sin(m*lon));
            q3 = q3 + m*pnm(n+1,m+1)*((*Global::Snm)(n+1,m+1)*cos(m*lon)-(*Global::Cnm)(n+1,m+1)*sin(m*lon));
        }
        dUdr     = dUdr     + q1*b1;
        dUdlatgc = dUdlatgc + q2*b2;
        dUdlon   = dUdlon   + q3*b3;
        q3 = 0; q2 = q3; q1 = q2;
    }
    
    // Body-fixed acceleration
    double r2xy = pow(r_bf(1), 2)+pow(r_bf(2),2);
    
    double ax = (1/d*dUdr-r_bf(3)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(1)-(1/r2xy*dUdlon)*r_bf(2);
    double ay = (1/d*dUdr-r_bf(3)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(2)+(1/r2xy*dUdlon)*r_bf(1);
    double az =  1/d*dUdr*r_bf(3)+sqrt(r2xy)/pow(d, 2)*dUdlatgc;

    Matrix a_bf(3);     a_bf(1)=ax;      a_bf(2) = ay;      a_bf(3) = az;
    // Inertial acceleration 
    return transpose(transpose(E)*transpose(a_bf));
}