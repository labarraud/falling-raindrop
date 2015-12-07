#ifndef FILE_TIME_SCHEME_HXX

#include "Matrix.cpp"


class VirtualOdeSystem
{
public:
  virtual ~VirtualOdeSystem();

  virtual void AddFunction(double alpha, const Matrix& rho, double t, Matrix& y) = 0;
};

#define FILE_TIME_SCHEME_HXX
#endif
