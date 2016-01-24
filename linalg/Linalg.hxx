#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <stdint.h>
#include <string>
#include <cstring>
#include <exception>

//! To display a variable (with its name)
#ifndef DISP
#define DISP(x) std::cout << #x ": " << x << std::endl
#endif

#include "Allocator.cxx"
#include "Vector.cxx"
#include "SparseVector.cxx"
#include "SparseMatrix.cxx"
#include "TinyVector.cxx"
#include "CoCg.cxx"
#include "CommonOutput.cxx"

#include "SolveMumps.cxx"

namespace linalg
{
  // typedef for points in R2 and R3
  typedef TinyVector<double, 2> R2;
  typedef TinyVector<double, 3> R3;
}
