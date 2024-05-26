// $Header$
//----------------------------------------------------------------------------------------
//                          matrix
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
/*
 * @file matrix.h
 * @class Matrix
 * @brief Class representing a matrix and its operations.
 *
 * @author Lorena Remacha Bordallo
*/

#ifndef _MATRIX_
#define _MATRIX_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;

class Matrix{
    public:
        int n_row, n_column;
        double **data;

        // Constructor parametrizado

        /**
         * @brief Parameterized constructor of vector.
         * @param v_size Size of the vector.
         */
        Matrix(const int v_size);
        /**
         * @brief Parameterized constructor of Matrix.
         * @param n_row Number of rows.
         * @param n_column Number of columns.
         */
        Matrix(const int n_row, const int n_column);

        // Operadores miembros
        /**
         * @brief Overloaded operator for accessing matrix elements.
         * @param n Index of the element (1D indexing).
         * @return Reference to the matrix element.
         */
        double& operator () (const int n);
        /**
         * @brief Overloaded operator for accessing matrix elements.
         * @param row Row index of the element.
         * @param column Column index of the element.
         * @return Reference to the matrix element.
         */
        double& operator () (const int row, const int column);
        /**
         * @brief Overloaded operator for scalar addition.
         * @param s Scalar value to add.
         * @return Resultant matrix.
         */
        Matrix& operator + (double s);
        /**
         * @brief Overloaded operator for matrix addition.
         * @param m Matrix to add.
         * @return Resultant matrix.
         */
        Matrix& operator + (Matrix &m);
        /**
         * @brief Overloaded operator for scalar subtraction.
         * @param s Scalar value to subtract.
         * @return Resultant matrix.
         */
        Matrix& operator - (double s);
        /**
         * @brief Overloaded operator for matrix subtraction.
         * @param m Matrix to subtract.
         * @return Resultant matrix.
         */
        Matrix& operator - (Matrix &m);
        /**
         * @brief Overloaded operator for scalar multiplication.
         * @param s Scalar value to multiply.
         * @return Resultant matrix.
         */
        Matrix& operator * (double s);
        /**
         * @brief Overloaded operator for matrix multiplication.
         * @param m Matrix to multiply.
         * @return Resultant matrix.
         */
        Matrix& operator * (Matrix &m);
        /**
         * @brief Overloaded operator for scalar division.
         * @param s Scalar value for division.
         * @return Resultant matrix.
         */
        Matrix& operator / (double s);

        // Operadores no miembros

        /**
         * @brief Overloaded operator for output stream.
         * @param o Output stream.
         * @param m Matrix to output.
         * @return Reference to the output stream.
         */
        friend ostream& operator << (ostream &o, Matrix &m);
};

// Sobrecarga de operadores
ostream& operator << (ostream &o, Matrix &m);

// Metodos
/**
 * @brief Function to generate a matrix of zeros.
 * @param n_row Number of rows.
 * @param n_column Number of columns.
 * @return Matrix of zeros.
 */
Matrix& zeros(const int n_row, const int n_column);
/**
 * @brief Function to generate an identity matrix.
 * @param size Size of the identity matrix.
 * @return Identity matrix.
 */
Matrix& eye(const int size); //devuelve matriz identidad de tamaño n
/**
 * @brief Function to transpose a matrix.
 * @param m Matrix to transpose.
 * @return Transposed matrix.
 */
Matrix& transpose(Matrix& m);
/**
 * @brief Function to compute the norm of a vector.
 * @param v Vector for which to compute the norm.
 * @return Norm of the vector.
 */
double norm(Matrix& v);
/**
 * @brief Function to compute the dot product of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Dot product of the two vectors.
 */
double dot(Matrix& v1, Matrix& v2);
/**
 * @brief Function to compute the cross product of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Cross product of the two vectors.
 */
Matrix& cross(Matrix& v1, Matrix& v2);
/**
 * @brief Function to compute the determinant of a square matrix.
 * @param v1 Matrix for which to compute the determinant.
 * @return Determinant of the matrix.
 */
double det(Matrix& v1);
/**
 * @brief Function to compute the inverse of a square matrix.
 * @param m Matrix to invert.
 * @return Inverted matrix.
 */
Matrix& inv(Matrix& m);
/**
 * @brief Function to extract a row from a matrix.
 * @param v1 Matrix from which to extract the row.
 * @param row Row index to extract.
 * @return Extracted row as a matrix.
 */
Matrix& extract_row(Matrix& v1, int row);
/**
 * @brief Function to extract a column from a matrix.
 * @param v1 Matrix from which to extract the column.
 * @param column Column index to extract.
 * @return Extracted column as a matrix.
 */
Matrix& extract_column(Matrix& v1, int column);
/**
 * @brief Function to assign a row to a matrix.
 * @param v1 Matrix to which to assign the row.
 * @param row Row index to assign.
 * @param v Row vector to assign.
 * @return Matrix with assigned row.
 */
Matrix& assign_row(Matrix& v1, int row, Matrix& v);
/**
 * @brief Function to assign a column to a matrix.
 * @param v1 Matrix to which to assign the column.
 * @param column Column index to assign.
 * @param v Column vector to assign.
 * @return Matrix with assigned column.
 */
Matrix& assign_column(Matrix& v1, int column, Matrix& v);
/**
 * @brief Function to convert a vector to an array.
 * @param v Vector to convert.
 * @return Pointer to the array.
 */
double* vectorToArray(Matrix& v) ;
/**
 * @brief Function to convert an array to a vector.
 * @param a Array to convert.
 * @param n Size of the array.
 * @return Converted vector.
 */
Matrix& arrayToVector(double a[], int n);


#endif