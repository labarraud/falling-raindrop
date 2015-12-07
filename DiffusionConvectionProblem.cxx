#ifndef DIFFUSION_CONVECTION_PROBLEM_CXX

#include "DiffusionConvectionProblem.hxx"
#include "TimeSchema.hxx"


class DiffusionConvectionProblem : public VirtualOdeSystem
{

protected:
	Vector<double> step_x, step_y;
	Matrix VX,VY;
	double D,delta_x,delta_t,delta_t;
	

public:
//AdvectionProblem();
	void Init(double delta_x, delta_t,double D);
	void ChampsCirculaire(double center, double intensite);
	void AddFunction(double alpha, const Matrix& rho, double t, Matrix& y);
	void WriteGnuPlot(Matrix M,const string& nom);
	
};


#define DIFFUSION_CONVECTION_PROBLEM_CXX
#endif
