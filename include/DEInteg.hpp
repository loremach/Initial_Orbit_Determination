// $Header$
//----------------------------------------------------------------------------------------
//                          DEInteg
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
 * @file DEInteg.hpp
 * @brief Header file for the ODE solver functions
 *
 * @details This header file contains the declarations for solving ordinary differential equations (ODEs) stepsize multistep method of Shampine & Gordon.
 *
 * @author Lorena Remacha Bordallo
*/

# include "matrix.h"

/**
 * @brief Carries out the ODE solution algorithm.
 *
 * This function implements a method to solve a system of ordinary differential equations (ODEs)
 * with adaptive step size and order control. It uses intermediate calculations to manage error
 * tolerance and step size adjustment dynamically.
 *
 * @param f User-supplied function that evaluates the right-hand side of the ODE.
 *          The function accepts input values T and Y[], evaluates the derivatives,
 *          and stores the result in YP[].
 * @param neqn Number of equations in the system.
 * @param y Array containing the current solution.
 * @param t Reference to the current value of the independent variable.
 * @param tout Desired value of T on output.
 * @param relerr Relative error tolerance.
 * @param abserr Absolute error tolerance.
 * @param iflag Reference to an integer indicating the status of integration.
 *              - 1: Normal input for initial call.
 *              - 2: Normal output, integration reached TOUT.
 *              - 3: Error tolerances too small, increased appropriately for continuing.
 *              - 4: More than 500 steps taken.
 *              - 5: Equations appear to be stiff.
 *              - 6: Invalid input parameters (fatal error).
 * @param yy Workspace array used to hold old solution data.
 * @param wt Error weight vector.
 * @param p Workspace array.
 * @param yp Workspace array used to hold values of the solution derivative.
 * @param ypout Workspace array used to hold values of the solution derivative.
 * @param phi Workspace array containing divided difference information about the polynomial interpolant to the solution.
 * @param alpha Workspace array.
 * @param beta Workspace array.
 * @param sig Workspace array.
 * @param v Workspace array.
 * @param w Workspace array.
 * @param g Workspace array.
 * @param phase1 Reference to a boolean indicating whether the program is in the first phase (always wants to increase the ODE method order).
 * @param psi Workspace array containing information about the polynomial interpolant to the solution.
 * @param x Reference to a "working copy" of T, the current value of the independent variable, which is adjusted as the code attempts to take a step.
 * @param h Reference to the current step size.
 * @param hold Reference to the last successful step size.
 * @param start Reference to a boolean that is TRUE on input for the first step. The program initializes data and sets START to FALSE.
 * @param told Reference to the previous value of T.
 * @param delsgn Reference to the sign (+1 or -1) of TOUT - T.
 * @param ns Reference to the number of steps taken with step size H.
 * @param nornd Reference to a boolean indicating whether rounding errors are negligible.
 * @param k Reference to the order of the current ODE method.
 * @param kold Reference to the order of the ODE method on the previous step.
 * @param isnold Reference to the previous value of ISN, the sign of IFLAG.
 */
void de ( Matrix (*f) (double, Matrix &), int neqn, double y[],
          double &t, double tout, double relerr, double abserr, int &iflag, double yy[],
          double wt[], double p[], double yp[], double ypout[], double phi[],
          double alpha[], double beta[], double sig[], double v[], double w[],
          double g[], bool &phase1, double psi[], double &x, double &h, double &hold,
          bool &start, double &told, double &delsgn, int &ns, bool &nornd, int &k, int &kold,
          int &isnold );

/**
 * @brief Returns the sign of an integer.
 *
 * @param i The integer whose sign is desired.
 * @return int The sign of the integer.
 */
int i4_sign ( int i );

/**
 * @brief Approximates the solution at XOUT by polynomial interpolation.
 *
 * @param x The point where the solution has been computed.
 * @param y The computed solution at X.
 * @param xout The point at which the solution is desired.
 * @param yout The solution at XOUT.
 * @param ypout The derivative of the solution at XOUT.
 * @param neqn The number of equations.
 * @param kold The order used for the last successful step.
 * @param phi Workspace array for interpolating polynomial information.
 * @param psi Workspace array for interpolating polynomial information.
 */
void intrp ( double x, double y[], double xout, double yout[], double ypout[],
             int neqn, int kold, double phi[], double psi[] );

/**
 * @brief Integrates a system of first order ODEs from T to TOUT.
 *
 * @param f User-supplied function that evaluates the right-hand sides of the ODE.
 * @param neqn The number of equations.
 * @param y The current solution vector.
 * @param t The current value of the independent variable.
 * @param tout The desired value of T on output.
 * @param relerr The relative error tolerance.
 * @param abserr The absolute error tolerance.
 * @param iflag Status indicator for the integration process.
 * @param work Workspace array.
 * @param iwork Workspace array for integer values.
 */
void ode ( Matrix (*f) (double, Matrix &), int neqn, double y[],
           double &t, double tout, double relerr, double abserr, int &iflag,
           double work[], int iwork[] );

/**
 * @brief Returns the sign of a double.
 *
 * @param x The number whose sign is desired.
 * @return double The sign of the number.
 */
double r8_sign ( double x );

/**
 * @brief Integrates the system of ODEs one step, from X to X+H.
 *
 * @param x The current value of the independent variable.
 * @param y The current solution vector.
 * @param f User-supplied function that evaluates the right-hand sides of the ODE.
 * @param neqn The number of equations.
 * @param h The current step size.
 * @param eps The local error tolerance.
 * @param wt The vector of error weights.
 * @param start Indicates whether this is the first step.
 * @param hold The last successful step size.
 * @param k The order of the current ODE method.
 * @param kold The order of the ODE method on the previous step.
 * @param crash Indicates whether no step can be taken.
 * @param phi Workspace array for divided difference information.
 * @param p Workspace array.
 * @param yp Workspace array for the solution derivative.
 * @param psi Workspace array for polynomial interpolant information.
 * @param alpha Workspace array.
 * @param beta Workspace array.
 * @param sig Workspace array.
 * @param v Workspace array.
 * @param w Workspace array.
 * @param g Workspace array.
 * @param phase1 Indicates whether the program is in the first phase.
 * @param ns The number of steps taken with step size H.
 * @param nornd Unused parameter.
 */
void step ( double &x, double y[], Matrix (*f) (double, Matrix &),
            int neqn, double &h, double &eps, double wt[], bool &start, double &hold,
            int &k, int &kold, bool &crash, double phi[], double p[], double yp[],
            double psi[], double alpha[], double beta[], double sig[], double v[],
            double w[], double g[], bool &phase1, int &ns, bool &nornd );

/**
 * @brief Prints the current YMDHMS date as a timestamp.
 */
void timestamp ( );

/**
 * @brief Computes the variational equations.
 *
 * @param f User-supplied function that evaluates the right-hand sides of the ODE.
 * @param neqn The number of equations.
 * @param y The current solution vector.
 * @param t The current value of the independent variable.
 * @param tout The desired value of T on output.
 * @param relerr The relative error tolerance.
 * @param abserr The absolute error tolerance.
 * @return Matrix The integrated result.
 */
Matrix DEInteg ( Matrix (*f) (double, Matrix &), int neqn, Matrix y,
           double t, double tout, double relerr, double abserr);
