#include <iostream>
#include <sstream>
#include <cmath>
#include "matrix.h"

using namespace matrix;

const int zeroSizeFlag(-1);

Matrix::Matrix(const int rows, const int columns): rows_(rows), columns_(columns)
{
    int size = rows_*columns_;
    if (size < 1) throw (zeroSizeFlag);

    data_ = new double[size];

    for (int i(0); i < size; i++) {
        data_[i] = 0;
    }
}

Matrix::Matrix(const int rows) : rows_(rows), columns_(1)
{
    int size = rows_*columns_;
    if (size < 1) throw (zeroSizeFlag);

    data_ = new double[size];

    for (int i(0); i < size; i++) {
        data_[i] = 0;
}

}Matrix::Matrix(const Matrix &rhs)
{
    if (rhs.getSize() > 0)
    {
        data_ = new double[rhs.getSize()];
        rows_ = rhs.rows_; columns_ = rhs.columns_;
        for (int i(0); i < getSize(); i++) {
            data_[i] = rhs.data_[i];
        }
    }
}

Matrix &Matrix::operator=(const Matrix &rhs)
{
    if (&rhs == this) return *this;

    delete[] data_;
    data_ = 0;
    rows_ = columns_ = 0;

    if (rhs.getSize() > 0) {
        data_ = new double[rhs.getSize()];
        rows_ = rhs.rows_;
        columns_ = rhs.columns_;
        for (int i(0); i < getSize(); i++) {
            data_[i] = rhs.data_[i];
            }
        }

    return *this;
}

double & Matrix::operator()(const int row, const int column) const
{
    if( (row < 1 || row > rows_) || (column < 1 || column > columns_) ) {
        std::stringstream errorOOB;
        errorOOB << "Matrix element out of bounds (" << row << ", " << column << ")";
        throw(errorOOB.str());
    }
    return data_[ (row - 1) * columns_ + (column - 1) ];
}

double & Matrix::operator[] (const int row) const
{
    if( (row < 1 || row > rows_) ) {
        std::stringstream errorOOB;
        errorOOB << "Matrix element out of bounds [" << row << "]";
        throw(errorOOB.str());
        //return 0;
    }
    return data_[ (row - 1) ];
}
Matrix Matrix::operator+(const Matrix &rhs) const
{
    if (rows_ != rhs.rows_ || columns_ != rhs.columns_) {
        throw("Matrices must have identical dimensions");
        return Matrix();
    }
    Matrix temp(rows_, columns_);
    for (int i(0); i < temp.getSize(); i++) {
        temp.data_[i] = data_[i] + rhs.data_[i];
    }
    return temp;
}

Matrix Matrix::operator-(const Matrix & rhs) const
{
    if (rows_ != rhs.rows_ || columns_ != rhs.columns_) {
        throw("Matrices must have identical dimensions");
        //return Matrix();
    }

    Matrix temp(rows_, columns_);
    for (int i(0); i < temp.getSize(); i++) {
        temp.data_[i] = data_[i] - rhs.data_[i];
    }
    return temp;
}

Matrix Matrix::operator*(const Matrix & rhs) const
{
    if (columns_ != rhs.rows_) {
        throw("Incompatible matrices under multiplication");
        //return Matrix();
        }

    Matrix temp(rows_, rhs.columns_);
    for (int n(1); !(n > temp.rows_); n++)
    {
        for (int m(1); !(m > temp.columns_); m++) // for each element of the temp Matrix, using indices starting at 1
        {
            double C(0); // initialise accumulator
            for (int i(1); !(i > columns_); i++)
            C += (*this)(n, i) * rhs(i, m);   // multiply M1(n,m) by M2(m,n)
            temp(n, m) = C;
        }
    }
    return temp;
}

Matrix Matrix::operator*(const double & rhs) const // scalar multiplication
{
    Matrix temp(rows_, columns_);
    for (int n(1); !(n > temp.rows_); n++)
    {
        for (int m(1); !(m > temp.columns_); m++) // for each element of the temp Matrix, using indices starting at 1
        {
            temp(n, m) = rhs * (*this)(n, m);
        }
    }
    return temp;
}

Matrix Matrix::transpose() const
{
    Matrix temp(columns_, rows_);
    for (int n(1); !(n > rows_); n++)
    {
        for (int m(1); !(m > columns_); m++)
        temp(m, n) = (*this)(n, m);
    }
    return temp;
}

double Matrix::det() const
{
    if (rows_ != columns_) {
        throw("Matrix must be square");
        //return 0;
        }

    if (!(rows_ == 2 || rows_ == 3)) {
        throw("Determinant function not valid for matrices larger than 3x3");
        //return 0;
        }

    double C(0), S_p, S_n;
    for (int i(1); !(i > columns_); i++) // move along first
    {
        S_p = S_n = (*this)(1, i); // set first row value
        for (int n(1); n < columns_; n++) // multiply down diagonals, Rule of Sarrus
        {
            if (i+n > columns_) S_p *= (*this)(1+n,(i+n)%columns_); // positive diagonal
            else S_p *= (*this)(1+n, i+n);

            if (i-n < 1) S_n *= (*this)(1+n,(columns_ - n + 1));
            else S_n *= (*this)(1+n, i-n);
        }
        C += S_p - S_n;
    }
    return C;
}

double Matrix::mod() const
{
    double accumulator(0);

    if (columns_ != 1) throw ("Matrix must be a vector");

    for (int n(1); n <= rows_; n++) {
        accumulator += (*this)[n] * (*this)[n];
    }

   return sqrt(accumulator);
}

std::ostream &matrix::operator<<(std::ostream &os, const Matrix &rhs)
{
    for (int n(0); n < rhs.rows_; n++)
    {
        for (int m(0); m < rhs.columns_; m++)
        {
            os << '\t';
            if (rhs.data_[n * rhs.columns_ + m] < 0) os << '\b'; // shift negatives numbers for better formatting
            os << rhs.data_[n * rhs.columns_ + m];
        }
        os << '\n';
    }
    os << std::endl;
    return os;
}
