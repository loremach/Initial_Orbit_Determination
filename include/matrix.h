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
        Matrix(const int v_size);
        Matrix(const int n_row, const int n_column);

        // Operadores miembros
        double& operator () (const int n);
        double& operator () (const int row, const int column);
        Matrix& operator + (double s);
        Matrix& operator + (Matrix &m);
        Matrix& operator - (double s);
        Matrix& operator - (Matrix &m);
        Matrix& operator * (double s);
        Matrix& operator * (Matrix &m);
        Matrix& operator / (double s);

        // Operadores no miembros
        friend ostream& operator << (ostream &o, Matrix &m);
};

// Sobrecarga de operadores
ostream& operator << (ostream &o, Matrix &m);

// Metodos
Matrix& zeros(const int n_row, const int n_column);
Matrix& eye(const int size); //devuelve matriz identidad de tamaño n
Matrix& transpose(Matrix& m);
double norm(Matrix& v);
double dot(Matrix& v1, Matrix& v2);
Matrix& cross(Matrix& v1, Matrix& v2);
double det(Matrix& v1);
Matrix& inv(Matrix& m);
Matrix& extract_row(Matrix& v1, int row);
Matrix& extract_column(Matrix& v1, int column);
Matrix& assign_row(Matrix& v1, int row, Matrix& v);
Matrix& assign_column(Matrix& v1, int column, Matrix& v);
double* vectorToArray(Matrix& v) ;
Matrix& arrayToVector(double a[], int n);


#endif