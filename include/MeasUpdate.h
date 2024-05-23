#ifndef _MEASUPDATE_
#define _MEASUPDATE_

#include <math.h>
#include "matrix.h"

void MeasUpdate(Matrix & K, Matrix & x, Matrix & P, double z, double g, double s, Matrix G, double n);

#endif