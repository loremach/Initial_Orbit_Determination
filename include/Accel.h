#ifndef _ACCEL_
#define _ACCEL_

#include <math.h>
#include "..\include\Matrix.h"
#include "..\include\IERS.h"
#include "..\include\timediff.h"
#include "..\include\PrecMatrix.h"
#include "..\include\NutMatrix.h"
#include "..\include\PoleMatrix.h"
#include "..\include\GHAMatrix.h"
#include "..\include\Mjday_TDB.h"
#include "..\include\JPL_Eph_DE430.h"
#include "..\include\AccelHarmonic.h"
#include "..\include\AccelPointMass.h"

Matrix Accel(double x, Matrix & Y);

#endif