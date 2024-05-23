//
// Created by lore0 on 21/05/2024.
//

#ifndef HGIBBS
#define HGIBBS

#include <math.h>
#include <string>
#include "..\include\Matrix.h"
#include "..\include\SAT_Const.h"
#include "..\include\unit.h"
#include "..\include\angl.h"

void hgibbs (Matrix r1, Matrix r2, Matrix r3, double Mjd1, double Mjd2, double Mjd3, Matrix & v2, double & theta, double &theta1, double & cop, string & error);
#endif
