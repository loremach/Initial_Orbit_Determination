//
// Created by lore0 on 23/05/2024.
//

#ifndef ANGLESG
#define ANGLESG

#include "Matrix.h"
#include "Geodetic.h"
#include "LTC.h"
#include "IERS.h"
#include "Global.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "SAT_const.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GhaMatrix.h"
#include "rpoly.h"
#include "gibbs.h"
#include "hgibbs.h"
#include "elements.h"

void anglesg (double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix Rs1, Matrix Rs2, Matrix Rs3, Matrix & r2, Matrix & v2);

#endif
