#include "..\include\EqnEquinox.h"

double EqnEquinox(double Mjd_TT)
{
    // Nutation in longitude and obliquity
    double dpsi;
    double deps;
    NutAngles(Mjd_TT, dpsi, deps);

    // Equation of the equinoxes
    return dpsi * cos(MeanObliquity(Mjd_TT));
}