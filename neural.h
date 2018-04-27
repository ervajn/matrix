#include <cmath>
#include "matrix.h"

#define ARG_CHECK(cond) if (!(cond)) throw std::invalid_argument(#cond " failed")

Matrix softmax(const Matrix& m)
{
  Matrix result(m.rows(), m.cols());
  std::transform(m.data().begin(), m.data().end(), result.data().begin(), [](double x){return std::exp(x);});
  const double s = std::accumulate(result.data().begin(), result.data().end(), 0.0);
  for (auto& x : result.data())
    x /= s;
  return result;
}
