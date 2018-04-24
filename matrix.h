#include <vector>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <numeric>
#include <cassert>

#define ARG_CHECK(cond) if (!(cond)) throw std::invalid_argument(#cond " failed")

class Matrix
{
public:
  explicit Matrix(std::size_t rows, std::size_t cols, double value=0)
    :cols_(cols)
    ,data_(cols * rows, value)
  {}

  explicit Matrix(std::size_t rows, std::size_t cols, const std::vector<double>& data)
    :cols_(cols)
    ,data_(data)
  {
    ARG_CHECK(data.size() == rows*cols);
  }

  static Matrix eye(std::size_t size)
  {
    Matrix m(size, size);
    for (std::size_t i=0; i<size; ++i)
      m(i,i) = 1;
    return m;
  }  

  const double& operator()(int row, int col) const { return data_[row*cols_ + col]; }
  double& operator()(int row, int col) { return data_[row*cols_ + col]; }
  const std::vector<double>& data() const { return data_; }
  
  std::size_t rows() const { return data_.size() / cols_; }
  std::size_t cols() const { return cols_; }
  std::tuple<std::size_t, std::size_t> dim() const { return std::make_tuple(rows(), cols()); }
  void resize(std::size_t newrows, std::size_t newcols)
  {
    Matrix mx(newrows, newcols);
    const std::size_t ncols = std::min(cols(), newcols);
    const std::size_t nrows = std::min(rows(), newrows);
    for (std::size_t r=0; r<nrows; ++r)
      for (std::size_t c=0; c<ncols; ++c)
        mx(r, c) = (*this)(r, c);
    swap(mx);    
  }   

  Matrix& iota(double start=0)
  {
    std::iota(data_.begin(), data_.end(), start);
    return *this;
  }

  Matrix trp() const
  {
    Matrix t(cols_, rows());
    for (std::size_t i=0; i<rows(); ++i)
      for (std::size_t j=0; j<cols(); ++j)
        t(j, i) = (*this)(i, j);
    return t;
  }

  Matrix& operator+=(const Matrix& rhs) { return combine(rhs, std::plus<double>()); } 
  Matrix& operator-=(const Matrix& rhs) { return combine(rhs, std::minus<double>()); } 

public:
  Matrix& operator*=(double value)
  {
    for (auto& x : data_)
      x *= value;
    return *this;    
  }
  
  Matrix& operator/=(double value)
  {
    for (auto& x : data_)
      x /= value;
    return *this;    
  } 

  double dot(const Matrix& rhs) const
  {
    return std::inner_product(data_.begin(), data_.end(), rhs.data_.begin(), 0);
  }

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
  
  void swap(Matrix& rhs)
  {
    std::swap(cols_, rhs.cols_);
    std::swap(data_, rhs.data_);
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
  
private:  
  std::size_t cols_;
  std::vector<double> data_;
};

std::ostream& operator<<(std::ostream& os, const Matrix& m) { return m.print(os); }

bool operator==(const Matrix& lhs, const Matrix& rhs) { return lhs.cols() == rhs.cols() && lhs.data() == rhs.data(); }
bool operator!=(const Matrix& lhs, const Matrix& rhs) { return !(lhs == rhs); }

Matrix operator+(const Matrix& lhs, const Matrix& rhs) { return Matrix(lhs) += rhs; }
Matrix operator-(const Matrix& lhs, const Matrix& rhs) { return Matrix(lhs) -= rhs; }
Matrix operator*(const Matrix& lhs, double value) { return Matrix(lhs) *= value; }
Matrix operator*(double value, const Matrix& rhs) { return Matrix(rhs) *= value; }
Matrix operator/(const Matrix& lhs, double value) { return Matrix(lhs) /= value; }

Matrix operator*(const Matrix& lhs, const Matrix& rhs)
{
  ARG_CHECK(lhs.cols() == rhs.rows());
  Matrix m(lhs.rows(), rhs.cols());
  for (std::size_t i=0; i<lhs.rows(); ++i)
    for (std::size_t j=0; j<rhs.cols(); ++j)
      for (std::size_t k=0; k<lhs.cols(); ++k)
        m(i, j) += lhs(i, k) * rhs(k, j);
  
  return m;
}  
