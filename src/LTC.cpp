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