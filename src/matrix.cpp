// $Source$
//----------------------------------------------------------------------------------------
//                          matrix
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
/*
 * @file matrix.cpp
 * @class Matrix
 * @brief Class representing a matrix and its operations.
 *
 * @author Lorena Remacha Bordallo
*/

#include "../include/matrix.h"

Matrix::Matrix(const int v_size)
{
    if (v_size <= 0)
    {
        cout << "Vector create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
    }
    this->n_row = 1;
    this->n_column = v_size;
    this->data = (double **)malloc(v_size * sizeof(double *));

    if (this->data == NULL)
    {
        cout << "Vector create: error in data\n";
        exit(EXIT_FAILURE);
    }

    this->data[0] = (double *)malloc(n_column * sizeof(double));

    for (int j = 0; j < v_size; j++)
    {
        this->data[0][j] = 0;
    }

}

Matrix::Matrix(const int n_row, const int n_column)
{
    if (n_row <= 0 || n_column <= 0)
    {
        cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
    }

    this->n_row = n_row;
    this->n_column = n_column;
    this->data = (double **)malloc(n_row * sizeof(double *));

    if (this->data == NULL)
    {
        cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n_row; i++)
    {
        this->data[i] = (double *)malloc(n_column * sizeof(double));

        for (int j = 0; j < n_column; j++)
        {
            this->data[i][j] = 0;
        }
    }
}

double &Matrix::operator()(const int row, const int column)
{
    if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column)
    {
        cout << "Matrix get: error in row/columns";
        exit(EXIT_FAILURE);
    }
    return this->data[row - 1][column - 1];
}

double &Matrix::operator()(const int n)
{
    if (n <= 0 || n > this->n_column)
    {
        cout << "Vector get: error in row/columns";
        exit(EXIT_FAILURE);
    }
    return this->data[0][n - 1];
}

Matrix &Matrix::operator+(double s)
{
    Matrix *result = new Matrix(this->n_row, this->n_column);

    for (int i = 0; i < this->n_row; i++)
        for (int j = 0; j < this->n_column; j++)
            (*result).data[i][j] = this->data[i][j] + s;

    return *result;
}

Matrix &Matrix::operator+(Matrix &m)
{
    if (m.n_row != this->n_row || m.n_column != this->n_column)
    {
        cout << "Matrix +: error in row/columns";
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(this->n_row, this->n_column);

    for (int i = 0; i < this->n_row; i++)
        for (int j = 0; j < this->n_column; j++)
            (*result).data[i][j] = this->data[i][j] + m.data[i][j];

    return *result;
}

Matrix &Matrix::operator-(double s)
{
    Matrix *result = new Matrix(this->n_row, this->n_column);

    for (int i = 0; i < this->n_row; i++)
        for (int j = 0; j < this->n_column; j++)
            (*result).data[i][j] = this->data[i][j] - s;

    return *result;
}

Matrix &Matrix::operator-(Matrix &m)
{
    if (m.n_row != this->n_row || m.n_column != this->n_column)
    {
        cout << "Matrix -: error in row/columns";
        exit(EXIT_FAILURE);
    }

    Matrix (*result) = new Matrix(this->n_row, this->n_column);

    for (int i = 0; i < this->n_row; i++)
        for (int j = 0; j < this->n_column; j++)
            (*result).data[i][j] = this->data[i][j] - m.data[i][j];

    return *result;
}

Matrix &Matrix::operator*(double s)
{
    Matrix *result = new Matrix(this->n_row, this->n_column);

    for (int i = 0; i < this->n_row; i++)
        for (int j = 0; j < this->n_column; j++)
            (*result).data[i][j] = this->data[i][j] * s;

    return *result;
}

Matrix& Matrix::operator*(Matrix &m)
{
    if (m.n_row != this->n_column)
    {
        cout << "Matrix *: error in row/columns";
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(this->n_row, m.n_column);

    for (int i = 0; i < this->n_row; i++)
    {
        for (int j = 0; j < m.n_column; j++)
        {
            (*result).data[i][j] = 0;
            for (int k = 0; k < this->n_column; k++)
            {
                (*result).data[i][j] = (*result).data[i][j] + this->data[i][k] * m.data[k][j];
            }
        }
    }
    return *result;
}

Matrix &Matrix::operator/(double s)
{

    if (s == 0.0)
    {
        cout << "Matrix /: invalid operation. Cannot be divided by zero";
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(this->n_row, this->n_column);

    for (int i = 0; i < this->n_row; i++)
        for (int j = 0; j < this->n_column; j++)
            (*result).data[i][j] = this->data[i][j] / s;

    return *result;
}

ostream &operator<<(ostream &o, Matrix &m)
{
    for (int i = 0; i < m.n_row; i++)
    {   
        for (int j = 0; j < m.n_column; j++)
        {
            o << m.data[i][j] << " ";
        }
        o << endl;
    }
    return o;
}

Matrix &zeros(const int n_row, const int n_column)
{
    if (n_row <= 0 || n_column <= 0)
    {
        cout << "Matrix zeros: error in row/columns";
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(n_row, n_column);
    return *result;
}

Matrix &eye(const int size)
{
    if (size <= 0)
    {
        cout << "Matrix eye: error in size";
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(size, size);

    for (int i = 0; i < size; i++)
        (*result).data[i][i] = 1;

    return *result;
}

Matrix &transpose(Matrix &m)
{
    Matrix *result = new Matrix(m.n_column, m.n_row);

    for (int i = 0; i < m.n_row; i++)
    {
        for (int j = 0; j < m.n_column; j++)
        {
            (*result).data[j][i] = m.data[i][j];
        }
    }

    return *result;
}

double norm(Matrix &v)
{
    double sum = 0.0;

    for (int i = 0; i < v.n_row; ++i)
    {
        for (int j = 0; j < v.n_column; ++j)
        {
            sum += v.data[i][j] * v.data[i][j];
        }
    }

    return sqrt(sum);
}

double dot(Matrix &v1, Matrix &v2)
{
    if (v1.n_row != 1 || v2.n_row != 1 || v1.n_column!=v2.n_column)
    {
        cout << "Matrix dot: invalid operation. Arrays with invalid number of components";
        exit(EXIT_FAILURE);
    }

    double res = 0.0;

    for (int i = 0; i < v1.n_column; i++)
    {
        res = res + v1.data[0][i] * v2.data[0][i];
    }
    return res;
}

Matrix &cross(Matrix &v1, Matrix &v2)
{
    if (v1.n_row != 1 || v2.n_row != 1 || v1.n_column != 3 || v2.n_column != 3)
    {
        cout << "Matrix cross: invalid operation. Arrays with invalid number of components";
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(3);

    (*result).data[0][0] = v1.data[0][1] * v2.data[0][2] - v1.data[0][2] * v2.data[0][1];
    (*result).data[0][1] = (-1) * (v1.data[0][0] * v2.data[0][2] - v1.data[0][2] * v2.data[0][0]);
    (*result).data[0][2] = v1.data[0][0] * v2.data[0][1] - v1.data[0][1] * v2.data[0][0];

    return *result;
}

double det(Matrix &v1)
{
    if ((v1.n_column != 1 || v1.n_row != 1) && (v1.n_column != 3 || v1.n_row != 3))
    {
        cout << "Matrix det: Matrix with invalid number of components";
        exit(EXIT_FAILURE);
    }

    if(v1.n_column == 3)
    {
        double det = v1.data[0][0] * (v1.data[1][1] * v1.data[2][2] - v1.data[1][2] * v1.data[2][1]) 
                    - v1.data[0][1] * (v1.data[1][0] * v1.data[2][2] - v1.data[1][2] * v1.data[2][0]) 
                    + v1.data[0][2] * (v1.data[1][0] * v1.data[2][1] - v1.data[1][1] * v1.data[2][0]);
        return det;
    }else if(v1.n_column == 1){
        double det = v1.data[0][0];
        return det;
    }else{
        cout << "Matrix det: Not implemented";
        exit(EXIT_FAILURE);
    }
}

Matrix &inv(Matrix &m)
{
    double determinant = det(m);
    if (determinant == 0.0)
    {
        printf("Matrix inv: The matrix has no inverse");
        exit(EXIT_FAILURE);
    }
    //Matrix *result;
    if(m.n_row==3){
        Matrix *result = new Matrix(3, 3);

        (*result).data[0][0] = (m.data[1][1] * m.data[2][2] - m.data[1][2] * m.data[2][1]) / determinant;
        (*result).data[0][1] = (m.data[0][2] * m.data[2][1] - m.data[0][1] * m.data[2][2]) / determinant;
        (*result).data[0][2] = (m.data[0][1] * m.data[1][2] - m.data[0][2] * m.data[1][1]) / determinant;

        (*result).data[1][0] = (m.data[1][2] * m.data[2][0] - m.data[1][0] * m.data[2][2]) / determinant;
        (*result).data[1][1] = (m.data[0][0] * m.data[2][2] - m.data[0][2] * m.data[2][0]) / determinant;
        (*result).data[1][2] = (m.data[0][2] * m.data[1][0] - m.data[0][0] * m.data[1][2]) / determinant;

        (*result).data[2][0] = (m.data[1][0] * m.data[2][1] - m.data[1][1] * m.data[2][0]) / determinant;
        (*result).data[2][1] = (m.data[0][1] * m.data[2][0] - m.data[0][0] * m.data[2][1]) / determinant;
        (*result).data[2][2] = (m.data[0][0] * m.data[1][1] - m.data[0][1] * m.data[1][0]) / determinant;

        return *result;
    }else if(m.n_row==1){
        Matrix *result = new Matrix(1, 1);
        (*result).data[0][0] = 1/m.data[0][0];
        return *result;
    }else{
        printf("Matrix inv: Not implemented");
        exit(EXIT_FAILURE);
    }
    //return *result;
}

Matrix &extract_row(Matrix &v1, int row)
{

    if (row > v1.n_row)
    {
        printf("Matrix extract_row: Invalid number of row");
        exit(EXIT_FAILURE);
    }
    Matrix *result = new Matrix(v1.n_column);

    for (int i = 0; i < v1.n_column; i++)
    {
        (*result).data[0][i] = v1.data[row-1][i];
    }

    return *result;
}

Matrix& extract_column(Matrix &v1, int column)
{

    if (column > v1.n_column)
    {
        printf("Matrix extract_column: Invalid number of column");
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(v1.n_row, 1);

    for (int i = 0; i < v1.n_row; i++)
    {
        (*result).data[i][0] = v1.data[i][column-1];
    }

    return *result;
}

Matrix& assign_row(Matrix &v1, int row, Matrix &v)
{
    if (row > v1.n_row || v1.n_column != v.n_column || v.n_row != 1)
    {
        printf("Matrix assign_row: Invalid number of row or invalid matrix");
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(v1.n_row, v1.n_column);

    for (int i = 0; i < v1.n_row; i++)
    {
        for (int j = 0; j < v1.n_column; j++)
        {
            if(i==row-1)
                (*result).data[i][j] = v.data[0][j];
            else
                (*result).data[i][j] = v1.data[i][j];
        }
    }
    return *result;
}

Matrix& assign_column(Matrix &v1, int column, Matrix &v)
{
    if (column > v1.n_column || v1.n_row != v.n_row || v.n_column != 1)
    {
        printf("Matrix assign_column: Invalid number of column or invalid matrix");
        exit(EXIT_FAILURE);
    }

    Matrix *result =new Matrix(v1.n_row, v1.n_column);

    for (int i = 0; i < v1.n_row; i++)
    {
        for (int j = 0; j < v1.n_column; j++)
        {
            if(j==column-1)
                (*result).data[i][j] = v.data[i][0];
            else
                (*result).data[i][j] = v1.data[i][j];
        }
    }
    return *result;
}

double* vectorToArray(Matrix& v) {
    double* a = new double[v.n_column];
    for (int i = 0; i < v.n_column; i++) {
        a[i] = v.data[0][i];
    }
    return a;
}

Matrix& arrayToVector(double a[], int n){
    Matrix *m = new Matrix(n);
    for(int i = 0; i<n; i++){
        (*m).data[0][i] = a[i];
    }
    return *m;
}