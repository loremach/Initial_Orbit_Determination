// $Header$
//----------------------------------------------------------------------------------------
//                          rpoly
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/05/20
//
/*
 * @file rpoly.h
 * @brief Jenkins-Traub real polynomial root finder implementation.
 *
 * @detail Translation of TOMS493 from FORTRAN to C. This implementation of Jenkins-Traub partially adapts
 * the original code to a C environment by restructuring many of the 'goto' controls to better fit
 * a block structured form. It also eliminates the global memory allocation in favor of local,
 * dynamic memory management. The calling conventions are slightly modified to return the number of roots found as the
 * function value.
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef RPOLY
#define RPOLY

#include "math.h"
#include "matrix.h"
/**
 * @brief Finds the roots of a polynomial using the Jenkins-Traub algorithm.
 *
 * This function implements the Jenkins-Traub algorithm (translated from FORTRAN to C)
 * to find the roots of a polynomial given its coefficients.
 *
 * @param op Pointer to the array of polynomial coefficients.
 * @param degree Degree of the polynomial.
 * @param zeror Pointer to the array where the real parts of the roots will be stored.
 * @param zeroi Pointer to the array where the imaginary parts of the roots will be stored.
 * @return The number of roots found.
 */
int rpoly(double *op, int degree, double *zeror, double *zeroi);
/**
 * @brief Calculates the zeros of the quadratic a*z^2 + b1*z + c.
 *
 * This function calculates the zeros of the quadratic a*z^2 + b1*z + c.
 *
 * @param a Coefficient of the quadratic term.
 * @param b1 Coefficient of the linear term.
 * @param c Constant term.
 * @param sr Pointer to store the real part of the smaller root.
 * @param si Pointer to store the imaginary part of the smaller root.
 * @param lr Pointer to store the real part of the larger root.
 * @param li Pointer to store the imaginary part of the larger root.
 */
void quad(double a,double b1,double c,double *sr,double *si,
          double *lr,double *li);
/**
 * @brief Computes up to L2 fixed shift k-polynomials, testing for convergence.
 *
 * This function performs variable-shift k-polynomial iteration for a quadratic factor, converging
 * only if the zeros are equimodular or nearly so.
 *
 * @param l2 The number of iterations to perform.
 * @param nz Pointer to an integer to store the number of zeros found.
 */
void fxshfr(int l2, int *nz);
/**
 * @brief Performs quadratic iteration for finding roots.
 *
 * This function performs quadratic iteration to find the roots of a polynomial.
 *
 * @param uu Pointer to the coefficient of the quadratic.
 * @param vv Pointer to the coefficient of the quadratic.
 * @param nz Pointer to an integer to store the number of zeros found.
 */
void quadit(double *uu,double *vv,int *nz);
/**
 * @brief Performs variable-shift H polynomial iteration for a real zero.
 *
 * This function performs variable-shift H polynomial iteration for a real zero.
 *
 * @param sss The starting iterate.
 * @param nz Pointer to an integer to store the number of zeros found.
 * @param iflag Pointer to an integer flag.
 */
void realit(double sss, int *nz, int *iflag);
/**
 * @brief Calculates scalar quantities used to compute the next k polynomial and new estimates of the quadratic coefficients.
 *
 * This function calculates scalar quantities used to compute the next k polynomial
 * and new estimates of the quadratic coefficients.
 *
 * @param type Pointer to an integer variable indicating how the calculations are normalized.
 */
void calcsc(int *type);
/**
 * @brief Computes the next k polynomials using scalars computed in calcsc.
 *
 * This function computes the next k polynomials using scalars computed in calcsc.
 *
 * @param type Pointer to an integer indicating the type of computation.
 */
void nextk(int *type);
/**
 * @brief Computes new estimates of the quadratic coefficients using the scalars computed in calcsc.
 *
 * This function computes new estimates of the quadratic coefficients using the scalars computed in calcsc.
 *
 * @param type Type of computation.
 * @param uu Pointer to the coefficient of the quadratic.
 * @param vv Pointer to the coefficient of the quadratic.
 */
void newest(int type,double *uu,double *vv);
/**
 * @brief Divides p by the quadratic 1,u,v placing the quotient in q and the remainder in a,b.
 *
 * This function divides p by the quadratic 1,u,v placing the quotient in q and the remainder in a,b.
 *
 * @param nn The degree of the polynomial p.
 * @param u The first coefficient of the quadratic divisor.
 * @param v The second coefficient of the quadratic divisor.
 * @param p Pointer to the array of polynomial coefficients.
 * @param q Pointer to the array where the quotient will be stored.
 * @param a Pointer to store the remainder.
 * @param b Pointer to store the remainder.
 */
void quadsd(int n,double *u,double *v,double *p,double *q,
            double *a,double *b);
/**
 * @brief Finds the roots of a polynomial.
 *
 * This function finds the roots of a polynomial given its coefficients using the
 * Jenkins-Traub algorithm implemented in the rpoly function.
 *
 * @param coef Matrix containing the coefficients of the polynomial.
 * @param real Matrix to store the real parts of the roots.
 * @param im Matrix to store the imaginary parts of the roots.
 */
void roots(Matrix coef, Matrix & real, Matrix & im);
#endif