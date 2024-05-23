#ifndef _IERS_
#define _IERS_

#include "matrix.h"
#include <math.h>
#include "SAT_Const.h"

void IERS(double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi,
          double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC, Matrix eop, 
          double Mjd_UTC, char interp = 'n');


#endif