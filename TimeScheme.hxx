#ifndef FILE_TIME_SCHEME_HXX



class VirtualTimeScheme
{
public:
  virtual ~VirtualTimeScheme();

  // fonctions qui retournent l'itere courant rho^n
  virtual Vector<Vector<double> >& GetIterate() = 0;
  virtual const Vector<Vector<double> >& GetIterate() const = 0;

  // fonction pour initialiser le schema en temps
  virtual void SetInitialCondition(double t0, double dt, Vector<Vector<double> >& rho0, VirtualOdeSystem& sys) = 0;

  // fonction principale qui avance le schema en temps
  virtual void Advance(int n, double tn, VirtualOdeSystem& sys) = 0;

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

  virtual Vector<Vector<double> >& GetIterate();
  virtual const Vector<Vector<double> >& GetIterate() const;

  // fonction pour initialiser le schema en temps
  virtual void SetInitialCondition(double t0, double dt_, Vector< Vector<double> >& rho0, VirtualOdeSystem& sys);

  // fonction principale qui avance le schema en temps
  virtual void Advance(int n, double tn, VirtualOdeSystem& sys);

};




#define FILE_TIME_SCHEME_HXX
#endif
