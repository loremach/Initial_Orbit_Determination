//
// Created by lore0 on 21/05/2024.
//

#ifndef VAREQN
#define VAREQN

#include <math.h>
#include "..\include\Matrix.h"
#include "..\include\global.h"
#include "..\include\IERS.h"
#include "..\include\timediff.h"
#include "..\include\PrecMatrix.h"
#include "..\include\NutMatrix.h"
#include "..\include\PoleMatrix.h"
#include "..\include\GHAMatrix.h"
#include "..\include\AccelHarmonic.h"
#include "..\include\G_AccelHarmonic.h"

Matrix VarEqn(double x, Matrix & yPhi);

#endif
