
#ifndef TIME_SCHEME_CXX
#define TIME_SCHEME_CXX

#include "TimeScheme.hxx"


//! Destructeur
VirtualTimeScheme::~VirtualTimeScheme()
{
}


/*************************
 * ExplicitEulerIterator *
 *************************/

//! default constructor
ExplicitEulerIterator::ExplicitEulerIterator()
{
  dt = 0;
}


//! retourne rho^n
Matrix& ExplicitEulerIterator::GetIterate()
{
  return rho;
}

//! retourne rho^n
const Matrix& ExplicitEulerIterator::GetIterate() const
{
  return rho;
}


//! fonction pour initialiser le schema en temps
void ExplicitEulerIterator::
SetInitialCondition(double t0, double dt_, Matrix& rho0, VirtualOdeSystem& sys)
{
  // on detruit rho0 quand on en a plus besoin
  dt = dt_;
  rho = rho0;
  rho0.Clear();
  rho_next = rho;
}


//! fonction principale qui avance le schema en temps
void ExplicitEulerIterator::Advance(int n, double tn, VirtualOdeSystem& sys)
{
  // Euler explicite : rho^n+1 = rho^n + dt f(t^n, rho^n)
  sys.AddFunction(dt, rho, tn, rho_next);
  rho = rho_next;
}


/********************************
 * LowStorageRungeKuttaIterator *
 ********************************/


//! default constructor
LowStorageRungeKuttaIterator::LowStorageRungeKuttaIterator()
{
  dt = 0;
}


//! libere la memoire utilisee
void LowStorageRungeKuttaIterator::Clear()
{
  rho.Clear();
  rho_next.Clear();
}


//! retourne l'itere rho^n
Matrix& LowStorageRungeKuttaIterator::GetIterate()
{
  return rho;
}


const Matrix& LowStorageRungeKuttaIterator::GetIterate() const
{
  return rho;
}


// fonction pour initialiser le schema en temps
void LowStorageRungeKuttaIterator
::SetInitialCondition(double t0, double dt_, Matrix& rho0, VirtualOdeSystem& sys)
{
  dt = dt_;
  rho = rho0;
  rho0.Clear();
  rho_next = rho;
}


// fonction principale qui avance le schema en temps
void LowStorageRungeKuttaIterator::Advance(int n, double tn, VirtualOdeSystem& sys)
{
  rho_next.Zero();
  sys.AddFunction(dt, rho, tn, rho_next);
  rho=rho + 1.496590219992291e-01*rho_next;
//  Add(1.496590219992291e-01, rho_next, rho);


 rho_next = -4.178904744998519e-01*rho_next;
  sys.AddFunction(dt, rho, tn+1.496590219992291e-01*dt, rho_next);
  rho=rho + 3.792103129996273e-01*rho_next;
//  Add(3.792103129996273e-01, rho_next, rho);


 rho_next= -1.192151694642677e+00*rho_next;
  sys.AddFunction(dt, rho, tn+3.704009573642048e-01*dt, rho_next);
  rho=rho + 8.229550293869817e-01*rho_next;
  //Add(8.229550293869817e-01, rho_next, rho);



  rho_next = -1.697784692471528e+00*rho_next;
  sys.AddFunction(dt, rho, tn+6.222557631344432e-01*dt, rho_next);
  rho=rho + 6.994504559491221e-01*rho_next;
 // Add(6.994504559491221e-01, rho_next, rho);

  rho_next = -1.514183444257156e+00*rho_next;
  sys.AddFunction(dt, rho, tn+9.582821306746903e-01*dt, rho_next);
  rho=rho + 1.530572479681520e-01*rho_next;
  //Add(1.530572479681520e-01, rho_next, rho);

  rho_next = rho;
}


#endif
