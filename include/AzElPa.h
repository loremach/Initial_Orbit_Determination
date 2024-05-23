#ifndef _AZELPA_
#define _AZELPA_

#include <math.h>
#include "matrix.h"
#include "SAT_const.h"

void AzElPa(Matrix s, double& Az, double& El, Matrix& dAds, Matrix& dEds);

#endif