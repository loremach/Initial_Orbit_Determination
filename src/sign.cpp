#include "..\include\sign_.h"

double sign_(double a, double b){
    double result;
    if (b>=0.0)
        result = abs(a);
    else
        result = - abs(a);

    return result;
}