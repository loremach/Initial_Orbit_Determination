#include "..\include\utils.h"

double custom_mod(double x, double y){
    return fmod(fmod(x, y) + y, y);
}