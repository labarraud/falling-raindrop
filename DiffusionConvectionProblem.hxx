#ifndef DIFFUSION_CONVECTION_PROBLEM_HXX

#include "TimeSchema.hxx"


class DiffusionConvectionProblem : public VirtualOdeSystem
{

protected:
	Vector<double> step_x, step_y;
	Velocity V;
	int Nx,Ny,Nt;
	double L,H,D,Delta_x,Delta_t,Delta_t;


public:
//DiffusionConvectionProblem();
	void Init(double delta_x, delta_t,double D);
	void ChampsCirculaire(double center, double intensite);
	void AddFunction(double alpha, const Matrix& rho, double t, Matrix& y);
	void WriteGnuPlot(Matrix M,const string& nom);

};


#define DIFFUSION_CONVECTION_PROBLEM_HXX
#endif
