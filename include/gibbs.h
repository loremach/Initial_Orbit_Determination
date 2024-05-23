#ifndef _GIBBS_
#define _GIBBS_

#include <math.h>
#include <string>
#include "..\include\Matrix.h"
#include "..\include\SAT_Const.h"
#include "..\include\unit.h"
#include "..\include\angl.h"

void gibbs(Matrix r1,Matrix r2,Matrix r3, Matrix & v2, double & theta, double & theta1, double & cop, string & error);

#endif