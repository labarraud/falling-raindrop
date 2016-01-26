#ifndef DIFFUSION_CONVECTION_PROBLEM_HXX

#include "TimeScheme.hxx"
#include "Velocity"

class DiffusionConvectionProblem : public VirtualOdeSystem
{

protected:
	Vector<double> step_x, step_y;
	Velocity* V;
	Particle n;
	int Nx,Ny,Nt;
	double L,H,D,Delta_x,Delta_y,Delta_t;


public:
	DiffusionConvectionProblem();
	virtual ~DiffusionConvectionProblem();
	void Init(int Nx, int Ny,double D);
	virtual void AddFunction(double alpha, const Vector<Vector<double> >& rho, double t, Vector<Vector<double> >& y) = 0;
	void WriteGnuPlot(const Vector<Vector<double> >& M,const string& nom);
	double GetX(int i) const;
	double GetY(int i) const;

};


#define DIFFUSION_CONVECTION_PROBLEM_HXX
#endif
