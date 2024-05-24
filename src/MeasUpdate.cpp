#include "..\include\MeasUpdate.h"

void MeasUpdate(Matrix & K, Matrix & x, Matrix & P, double z, double g, double s, Matrix G, double n){
    double m=1;

    Matrix Inv_W = zeros(m,m);
    Inv_W(1, 1) = s*s;    // Inverse weight (measurement covariance)

    // Kalman gain
    K = transpose(P*transpose(G)*inv(Inv_W+G*P*transpose(G)));

    // State update
    x = x + K*(z-g);

    // Covariance update
    P = (eye(n)-transpose(K)*G)*P;
}