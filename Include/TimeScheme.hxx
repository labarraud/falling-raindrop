#ifndef FILE_TIME_SCHEME_HXX

#include "Matrix.hxx"
#include "DiffusionConvectionProblem.hxx"
#include <string>
#include <cstring>
#include <vector>



class VirtualTimeScheme
{
public:
  virtual ~VirtualTimeScheme();

  // fonctions qui retournent l'itere courant rho^n
  virtual Matrix& GetIterate() = 0;
  virtual const Matrix& GetIterate() const = 0;

  // fonction pour initialiser le schema en temps
  virtual void SetInitialCondition(double t0, double dt, Matrix& rho0, VirtualOdeSystem& sys) = 0;

  // fonction principale qui avance le schema en temps
  virtual void Advance(int n, double tn, VirtualOdeSystem& sys, const vector<precision>& sec_membre) = 0;

};




//! Schema de Runge-Kutta avec low-storage
class LowStorageRungeKuttaIterator : public VirtualTimeScheme
{
private:
  // on stocke le pas de temps
  double dt;
  // pour ce schema, on a besoin de ne stocker que deux vecteurs
  Matrix rho, rho_next;

public:
  LowStorageRungeKuttaIterator();

  void Clear();

  virtual Matrix& GetIterate();
  virtual const Matrix& GetIterate() const;

  // fonction pour initialiser le schema en temps
  virtual void SetInitialCondition(double t0, double dt_, Matrix& rho0, VirtualOdeSystem& sys);

  // fonction principale qui avance le schema en temps
  virtual void Advance(int n, double tn, VirtualOdeSystem& sys, const vector<precision>& sec_membre);

};


//! Exemple de classe derivee de VirtualTimeScheme : Euler explicite
class ExplicitEulerIterator : public VirtualTimeScheme
{
private:
  // on stocke le pas de temps
  double dt;
  // pour le schema d'Euler, on a besoin de ne stocker que deux vecteurs
  Matrix rho, rho_next;

public:
  ExplicitEulerIterator();

  virtual Matrix& GetIterate();
  virtual const Matrix& GetIterate() const;

  // fonction pour initialiser le schema en temps
  virtual void SetInitialCondition(double t0, double dt_, Matrix& rho0, VirtualOdeSystem& sys);

  // fonction principale qui avance le schema en temps
  virtual void Advance(int n, double tn, VirtualOdeSystem& sys, const vector<precision>& sec_membre);

};

void error_orderxy_circle(precision mindxy,precision hdxy,precision maxdxy
		,precision cfl ,precision tmaxdemi,precision omega,DiffusionConvectionProblem& ode,
		VirtualTimeScheme& time, const string& fileout);


//! Schema de Runge-Kutta d'ordre 2
class RK2Iterator : public VirtualTimeScheme
{
private:
  // on stocke le pas de temps
  double dt;
  // pour ce schema, on a besoin de ne stocker que deux vecteurs
  Matrix rho, rho_next;

public:
  RK2Iterator();

  void Clear();

  virtual Matrix& GetIterate();
  virtual const Matrix& GetIterate() const;

  // fonction pour initialiser le schema en temps
  virtual void SetInitialCondition(double t0, double dt_, Matrix& rho0, VirtualOdeSystem& sys);

  // fonction principale qui avance le schema en temps
  virtual void Advance(int n, double tn, VirtualOdeSystem& sys, const vector<precision>& sec_membre);

};



#define FILE_TIME_SCHEME_HXX
#endif
