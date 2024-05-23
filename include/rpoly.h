//
// Created by lore0 on 23/05/2024.
//

#ifndef RPOLY
#define RPOLY

#include "math.h"
#include "matrix.h"

void quad(double a,double b1,double c,double *sr,double *si,
          double *lr,double *li);
void fxshfr(int l2, int *nz);
void quadit(double *uu,double *vv,int *nz);
void realit(double sss, int *nz, int *iflag);
void calcsc(int *type);
void nextk(int *type);
void newest(int type,double *uu,double *vv);
void quadsd(int n,double *u,double *v,double *p,double *q,
            double *a,double *b);
void roots(Matrix coef, Matrix & real, Matrix & im);
#endif