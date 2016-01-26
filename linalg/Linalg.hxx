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


