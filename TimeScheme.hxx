#ifndef FILE_TIME_SCHEME_HXX

#include "Matrix.cpp"


class VirtualOdeSystem
{
public:
  virtual ~VirtualOdeSystem();

  virtual void AddFunction(double alpha, const Matrix& rho, double t, Matrix& y) = 0;
};


//! Schema de Runge-Kutta avec low-storage
class LowStorageRungeKuttaIterator : public VirtualTimeScheme
{
private:
  // on stocke le pas de temps
  double dt;
  // pour ce schema, on a besoin de ne stocker que deux vecteurs
  Vector< Vector<double> > rho, rho_next;

public:
  LowStorageRungeKuttaIterator();

  void Clear();

  virtual Vector<double>& GetIterate();
  virtual const Vector<double>& GetIterate() const;

  // fonction pour initialiser le schema en temps
  virtual void SetInitialCondition(double t0, double dt_, Vector<double>& rho0, VirtualOdeSystem& sys);

  // fonction principale qui avance le schema en temps
  virtual void Advance(int n, double tn, VirtualOdeSystem& sys);

};




#define FILE_TIME_SCHEME_HXX
#endif
