// $Source$
//----------------------------------------------------------------------------------------
//                          Legendre
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/14
//
/*
 * @file Legendre.cpp
 * @brief Computation of Legendre polynomials and their derivatives
 *
 * @details This file contains the implementation of functions for computing Legendre
 * polynomials and their derivatives.
 *
 * @author Lorena Remacha Bordallo
*/

#include "..\include\Legendre.h"

void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm){
    int aux;
    if (n>m)
        aux=n;
    else
        aux=m;

    pnm = zeros(n+1,aux+1);
    dpnm = zeros(n+1,aux+1);

    pnm(1,1)=1;
    dpnm(1,1)=0;
    pnm(2,2)=sqrt(3.0)*cos(fi);
    dpnm(2,2)=-sqrt(3.0)*sin(fi);

    // diagonal coefficients
    for (int i = 2; i <= n; i++){    
        pnm(i+1,i+1)= sqrt((2.0*i+1)/(2.0*i))*cos(fi)*pnm(i,i);
    }
    for (int i = 2; i <= n; i++){
        dpnm(i+1,i+1)= sqrt((2.0*i+1)/(2.0*i))*((cos(fi)*dpnm(i,i))-(sin(fi)*pnm(i,i)));
    }
    
    // horizontal first step coefficients
    for (int i = 1; i <= n; i++){
        pnm(i+1,i)= sqrt(2.0*i+1)*sin(fi)*pnm(i,i);
    }
    for (int i = 1; i <= n; i++){
        dpnm(i+1,i)= sqrt(2.0*i+1)*((cos(fi)*pnm(i,i))+(sin(fi)*dpnm(i,i)));
    }
    
    // horizontal second step coefficients
    int j=0;
    int k=2;
    while(true){
        for (int i = k; i <= n; i++){      
            pnm(i+1,j+1)=sqrt((2.0*i+1)/((i-j)*(i+j)))*((sqrt(2.0*i-1)*sin(fi)*pnm(i,j+1))-(sqrt(((i+j-1)*(i-j-1))/(2.0*i-3))*pnm(i-1,j+1)));
        }
        j = j+1;
        k = k+1;
        if (j>m){
            break;
        }
    }
    j = 0;
    k = 2;
    while(true){
        for (int i = k; i <= n; i++){       
            dpnm(i+1,j+1)=sqrt((2.0*i+1)/((i-j)*(i+j)))*((sqrt(2.0*i-1)*sin(fi)*dpnm(i,j+1))+(sqrt(2.0*i-1)*cos(fi)*pnm(i,j+1))-(sqrt(((i+j-1)*(i-j-1))/(2.0*i-3))*dpnm(i-1,j+1)));
        }
        j = j+1;
        k = k+1;
        if (j>m){
            break;
        }
    }
}