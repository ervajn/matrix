// g++ --std=c++11 main.cpp -Wall
#include <iostream>
#include <cassert>
#include <sstream>
#include "matrix.h"

static void testMisc()
{
  Matrix m1(10, 10, 1);
  Matrix m2(10, 10, 2);

  m1 += m2;
  m1 -= m2;

  assert(m1.dot(m2) == 200);

  Matrix a(3, 2, 1), b(2, 1, 2); 
  assert(a * b == Matrix(3, 1, 4));

  Matrix m3(10,10);
  m3.iota();
  assert(m3(0,0) == 0);
  assert(m3(9,9) == 99);  
  
  Matrix m(2, 3);
  m.iota(1);
  const Matrix mt(3, 2, std::vector<double>({1, 4, 2, 5, 3, 6}));
  assert(m.trp() == mt);
  
  assert(2 * Matrix::eye(5) * 2 / 4 == Matrix::eye(5));  

  Matrix m5(5, 5, 9);
  m5.resize(7, 7);
  assert(m5(6, 6) == 0);  
  m5.resize(4, 3);
  assert(m5 == Matrix(4, 3, 9));
}

static void testConvolute()
{
  const Matrix matrix(5, 6, 1);
  const Matrix filter(2, 2, 2);
  const Matrix result = matrix.convolute(filter);
  assert(result == Matrix(4, 5, 8));  
  // std::cout << "Input:\n" << matrix << std::endl;
  // std::cout << "Filter:\n" << filter << std::endl;
  // std::cout << "Result:\n" << result << std::endl;  
}

static void testIO()
{
  Matrix m(3, 4);
  m.iota();
  std::cout << m << std::endl;  
}


int main(int argc, char* argv[])
{
  testMisc();
  testConvolute();
  testIO();

  return 0;
}
