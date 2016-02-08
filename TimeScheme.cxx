
#ifndef TIME_SCHEME_CXX
#define TIME_SCHEME_CXX

#include "TimeScheme.hxx"

//! Destructeur
VirtualTimeScheme::~VirtualTimeScheme()
{
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
Vector<Vector<double> >& LowStorageRungeKuttaIterator::GetIterate()
{
  return rho;
}


const Vector<Vector<double> >& LowStorageRungeKuttaIterator::GetIterate() const
{
  return rho;
}


// fonction pour initialiser le schema en temps
void LowStorageRungeKuttaIterator
::SetInitialCondition(double t0, double dt_, Vector<Vector<double> >& rho0, VirtualOdeSystem& sys)
{
  dt = dt_;
  rho = rho0;
  rho0.Clear();
  rho_next = rho;
}


// fonction principale qui avance le schema en temps
void LowStorageRungeKuttaIterator::Advance(int n, double tn, VirtualOdeSystem& sys)
{
  for(int i(0); i < rho_next.GetM(); ++i) {
	  rho_next(i).Zero();
  }
  sys.AddFunction(dt, rho, tn, rho_next);
  Add(1.496590219992291e-01, rho_next, rho);

  for(int i(0); i < rho_next.GetM(); ++i) {
	  rho_next(i) *= -4.178904744998519e-01;
  }
  sys.AddFunction(dt, rho, tn+1.496590219992291e-01*dt, rho_next);
  Add(3.792103129996273e-01, rho_next, rho);

  for(int i(0); i < rho_next.GetM(); ++i) {
	  rho_next(i) *= -1.192151694642677e+00;
  }
  sys.AddFunction(dt, rho, tn+3.704009573642048e-01*dt, rho_next);
  Add(8.229550293869817e-01, rho_next, rho);

  for(int i(0); i < rho_next.GetM(); ++i) {
	  rho_next(i) *= -1.697784692471528e+00;
  }
  sys.AddFunction(dt, rho, tn+6.222557631344432e-01*dt, rho_next);
  Add(6.994504559491221e-01, rho_next, rho);

  for(int i(0); i < rho_next.GetM(); ++i) {
	  rho_next(i) *= -1.514183444257156e+00;
  }
  sys.AddFunction(dt, rho, tn+9.582821306746903e-01*dt, rho_next);
  Add(1.530572479681520e-01, rho_next, rho);

  rho_next = rho;
}


#endif
