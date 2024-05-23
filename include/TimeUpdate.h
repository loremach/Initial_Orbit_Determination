#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include <math.h>
#include "..\include\Matrix.h"

Matrix TimeUpdate(Matrix P, Matrix Phi, double Qdt=0.0);

#endif