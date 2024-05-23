#include "..\include\PoleMatrix.h"

Matrix PoleMatrix(double xp, double yp){
    Matrix rx = R_x(-yp);
    Matrix ry = R_y(-xp);
    
    return ry * rx;
}