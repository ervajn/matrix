#ifndef _INC_MATRIX_180427
#define _INC_MATRIX_180427

#include <vector>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <numeric>
#include <cassert>
#include <initializer_list>

#define ARG_CHECK(cond) if (!(cond)) throw std::invalid_argument(#cond " failed")

class Matrix;
std::ostream& operator<<(std::ostream& os, const Matrix& m);

class Matrix
{
public:
  using value_type = double;
  using data_type = std::vector<value_type>;
  using dim_type = std::tuple<std::size_t, std::size_t>;

  explicit Matrix(std::size_t rows, std::size_t cols, const data_type& data)
    :cols_(cols)
    ,data_(data)
  {
    ARG_CHECK(data.size() == rows*cols);
  }

  explicit Matrix(std::size_t rows=0, std::size_t cols=1, double value=0)
    :Matrix(rows, cols, data_type(cols * rows, value))
  {}

  explicit Matrix(std::size_t rows, std::size_t cols, const std::initializer_list<value_type>& l)
    :Matrix(rows, cols, data_type(l))
  {}

  static Matrix eye(std::size_t size)
  {
    Matrix m(size, size);
    for (std::size_t i=0; i<size; ++i)
      m(i,i) = 1;
    return m;
  }  

  const double& operator()(int row, int col) const { return data_[row*cols_ + col]; }
  double& operator()(int row, int col) { return data_[row*cols_ + col]; }
  const data_type& data() const { return data_; }
  data_type& data() { return data_; }
  
  std::size_t rows() const { return data_.size() / cols_; }
  std::size_t cols() const { return cols_; }
  dim_type dim() const { return std::make_tuple(rows(), cols()); }
  void resize(const dim_type& dim) { resize(std::get<0>(dim), std::get<1>(dim)); }  
  void resize(std::size_t newrows, std::size_t newcols)
  {
    if (newrows != rows() || newcols != cols())
    {
      Matrix mx(newrows, newcols);
      const std::size_t ncols = std::min(cols(), newcols);
      const std::size_t nrows = std::min(rows(), newrows);
      for (std::size_t r=0; r<nrows; ++r)
        for (std::size_t c=0; c<ncols; ++c)
          mx(r, c) = (*this)(r, c);
      swap(mx);
    }
  }x

  Matrix& iota(double start=0) { std::iota(data_.begin(), data_.end(), start); return *this; }

  Matrix& trp()
  {
    Matrix t(cols_, rows());
    for (std::size_t i=0; i<rows(); ++i)
      for (std::size_t j=0; j<cols(); ++j)
        t(j, i) = (*this)(i, j);
    swap(t);
    return *this;
  }

  Matrix& operator+=(const Matrix& rhs) { return combine(rhs, std::plus<double>()); } 
  Matrix& operator-=(const Matrix& rhs) { return combine(rhs, std::minus<double>()); } 
  Matrix& operator*=(double value) { for (auto& x : data_) x *= value; return *this; }
  Matrix& operator/=(double value) { for (auto& x : data_) x /= value; return *this; }
  double dot(const Matrix& rhs) const { return std::inner_product(data_.begin(), data_.end(), rhs.data_.begin(), 0.0); }

  void swap(Matrix& rhs) { std::swap(cols_, rhs.cols_); std::swap(data_, rhs.data_); }

  std::ostream& print(std::ostream& os) const
  {
    for (std::size_t i=0; i<rows(); ++i)
    {
      os << (i?"\n":"");
      for (std::size_t j=0; j<cols(); ++j)
        os << (j?" ":"") << std::setw(5) << (*this)(i,j);
    }
    return os;
  }

  Matrix convolute(const Matrix& filter) const
  {
    ARG_CHECK(filter.rows() <= rows());
    ARG_CHECK(filter.cols() <= cols());
    Matrix result(rows() - filter.rows() + 1, cols() - filter.cols() + 1);
    for (std::size_t r=0; r<result.rows(); ++r)
      for (std::size_t c=0; c<result.cols(); ++c)
        result(r, c) = applyFilter(filter, r, c);
    return result;
  }

  Matrix& cross(const Matrix& lhs, const Matrix& rhs)
  {
    ARG_CHECK(rows() == lhs.rows());
    ARG_CHECK(cols() == rhs.cols());
    ARG_CHECK(lhs.cols() == rhs.rows());
    for (std::size_t i=0; i<lhs.rows(); ++i)
      for (std::size_t j=0; j<rhs.cols(); ++j)
        for (std::size_t k=0; k<lhs.cols(); ++k)
          (*this)(i, j) += lhs(i, k) * rhs(k, j);
    return *this;
  }

  // Ax = B
  void solve(Matrix& A, Matrix& B)
  {
    ARG_CHECK(A.rows() == B.rows());
    const std::size_t n = A.rows();
    for (std::size_t k=0; k<n-1; ++k)
    {
      const std::size_t v = A.pivotRow(k);
      A.swapRows(k, v);
      B.swapRows(k, v);
      for (std::size_t i=k+1; i<n; ++i)
      {
        const value_type mik = A(i, k) / A(k, k);
        for (std::size_t j = k+1; j<n; ++j)
        {
          A(i, j) -= mik * A(k, j);
        }
        B(i, 0) -= mik * B(k, 0);
      }
    }
    resize(n, 1);
    for (int i=n-1; i>=0; --i)
    {
      double sum = B(i, 0);
      for (std::size_t j=i+1; j<n; ++j)
        sum -= A(i, j) * (*this)(j, 0);
      (*this)(i, 0) = sum / A(i, i);
    }
  }

private:
  double applyFilter(const Matrix& filter, std::size_t row, std::size_t col) const
  {
    assert(row + filter.rows() <= rows());
    assert(col + filter.cols() <= cols());
    double result = 0;
    for (std::size_t r=0; r<filter.rows(); ++r)
      for (std::size_t c=0; c<filter.cols(); ++c)
        result += (*this)(row + r, col + c) * filter(r, c);
    return result;
  }

  template <typename OP>
  Matrix& combine(const Matrix& rhs, OP op)
  {
    ARG_CHECK(dim() == rhs.dim());
    std::transform(data_.begin(), data_.end(), rhs.data_.begin(), data_.begin(), op);
    return *this;
  }

  void swapRows(std::size_t row1, std::size_t row2)
  {
    if (row1 == row2)
      return;
    const auto b1 = data_.begin() + row1 * cols();
    const auto e1 = b1 + cols();
    const auto b2 = data_.begin() + row2 * cols();
    std::swap_ranges(b1, e1, b2);
  }

  std::size_t pivotRow(std::size_t k) const
  {
    std::size_t maxRow = 0;
    value_type maxValue = 0.0;
    for (std::size_t r = k; r < rows(); ++r)
    {
      const value_type tmp = std::abs((*this)(r, k));
      if (tmp > maxValue)
      {
        maxValue = tmp;
        maxRow = r;
      }
    }
    return maxRow;
  }

private:
  std::size_t cols_;
  std::vector<double> data_;
};

inline std::ostream& operator<<(std::ostream& os, const Matrix& m) { return m.print(os); }

inline bool operator==(const Matrix& lhs, const Matrix& rhs) { return lhs.cols() == rhs.cols() && lhs.data() == rhs.data(); }
inline bool operator!=(const Matrix& lhs, const Matrix& rhs) { return !(lhs == rhs); }

inline Matrix operator+(const Matrix& lhs, const Matrix& rhs) { return Matrix(lhs) += rhs; }
inline Matrix operator-(const Matrix& lhs, const Matrix& rhs) { return Matrix(lhs) -= rhs; }
inline Matrix operator*(const Matrix& lhs, double value) { return Matrix(lhs) *= value; }
inline Matrix operator*(double value, const Matrix& rhs) { return Matrix(rhs) *= value; }
inline Matrix operator/(const Matrix& lhs, double value) { return Matrix(lhs) /= value; }
inline Matrix operator*(const Matrix& lhs, const Matrix& rhs)  { return Matrix(lhs.rows(), rhs.cols()).cross(lhs, rhs); }

#endif
