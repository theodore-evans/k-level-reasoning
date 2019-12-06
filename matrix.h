#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

namespace matrix
{
    class Matrix
    {
        friend std::ostream & operator<<(std::ostream & os, const Matrix &rhs);

        private:
            int rows_, columns_;
            double *data_;

        public:
            Matrix() : rows_(0), columns_(0), data_(0) {}
            Matrix(const int rows, const int columns);// prototype parameterised constructor
            Matrix(const int rows);

            Matrix & operator=(const Matrix &rhs);
            Matrix(const Matrix &rhs);

            ~Matrix() {delete[] data_;}

            int rows() const {return rows_;}
            int columns() const {return columns_;}
            int getSize() const {return rows_ * columns_;}
            void print() const;

            double &operator()(const int row, const int column) const;
            double &operator[] (const int row) const;

            Matrix operator+(const Matrix &rhs) const;
            Matrix operator-(const Matrix &rhs) const;
            Matrix operator*(const Matrix &rhs) const;
            Matrix operator*(const double &rhs) const;
            Matrix transpose() const;
            double det() const;
            double mod() const;
    };

    std::ostream & operator<<(std::ostream & os, const Matrix &rhs);
}

#endif // MATRIX_H_INCLUDED

