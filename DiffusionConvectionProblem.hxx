#ifndef DIFFUSION_CONVECTION_PROBLEM_HXX

#include "velocity.hxx"


//! Classe abstraite pour la definition d'une EDO
/*!
  une EDO s'ecrit rho' = f(t, rho)  
 */
class VirtualOdeSystem
{
public:
  virtual ~VirtualOdeSystem();

  // fonction pour calculer y = y + alpha f(t, rho)
  virtual void AddFunction(double alpha, const Vector<Vector<double> >& rho, double t, Vector<Vector<double> >& y) = 0;
  
};

class DiffusionConvectionProblem : public VirtualOdeSystem
{

protected:
	Vector<double> step_x, step_y;
	Velocity velocity;
	Particle particule;
	int Nx,Ny,Nt;
	double L,H,tfinal,D,Delta_x,Delta_y,Delta_t;


public:
	DiffusionConvectionProblem(int Nx,int Ny,int Nt,double L,double H,double tfinal,Velocity& V,Particle& n);
	virtual ~DiffusionConvectionProblem();
	//void Init(int Nx, int Ny,double D);
	virtual void AddFunction(double alpha, const Vector<Vector<double> >& rho, double t, Vector<Vector<double> >& y) = 0;
	void WriteGnuPlot(const Vector<Vector<double> >& M,const string& nom);
	double GetX(int i) const;
	double GetY(int i) const;
	Particle & GetP();

};


#define DIFFUSION_CONVECTION_PROBLEM_HXX
#endif
