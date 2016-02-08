#ifndef DIFFUSION_CONVECTION_PROBLEM_CXX

#include "DiffusionConvectionProblem.hxx"
#include "TimeScheme.hxx"


//! Destructeur
VirtualOdeSystem::~VirtualOdeSystem()
{
}




/*	Vector<double> step_x, step_y;
	Velocity V;
	Particle particule;
	int Nx,Ny,Nt;
	double L,H,D,Delta_x,Delta_y,Delta_t;*/

	DiffusionConvectionProblem::DiffusionConvectionProblem(int Nx,int Ny,int Nt,double L,double H,double tfinal,Velocity& V,Particle& n)
		:	velocity(V), particule(n)
	{

		this->Nx=Nx;
		this->Ny=Ny;
		this->L=L;
    	this->H=H;
    	this->D=0;
    	this->tfinal=H;
	    this->Delta_x=L/Nx;
		this->Delta_y=H/Ny;
		this->Delta_t=tfinal/Nt;


	}

	DiffusionConvectionProblem::~DiffusionConvectionProblem()
	{

	}

	/*void DiffusionConvectionProblem::Init(int Nx, int Ny,double D)
	{
		step_x.Reallocate(Nx+1);
		step_y.Reallocate(Ny+1);
		Delta_x = D/Nx;
		Delta_y = D/Ny;
		L = H = this.D = D;
		if(V != NULL) {
			delete V;
		}
		V = new Velocity(Nx, Ny, L, H);
		for(int i(0); i < Nx+1; ++i) {
			step_x(i) = i*Delta_x;
		}
		for(int i(0); i < Ny+1; ++i) {
			step_y(i) = i*Delta_y;
		}
	}
*/

	void DiffusionConvectionProblem::WriteGnuPlot(const Vector<Vector<double> >& M,const string& nom)
	{
		ofstream file_out(nom.data());
		file_out.precision(15);
		for(int i(0), j; i < M.GetM(); i++) {
			for(j = 0; j < M(0).GetM(); ++j) {
				file_out << step_x(j) << " " << step_y(i) << " " << M(i)(j) << '\n';
			}
			file_out  << '\n';
		}

		file_out.close();
	}

	double DiffusionConvectionProblem::GetX(int i) const
	{
		return step_x(i);
	}

	double DiffusionConvectionProblem::GetY(int i) const
	{
		return step_y(i);
	}
	
	Particle & DiffusionConvectionProblem::GetP()
	{
		return particule;
		
	}

#define DIFFUSION_CONVECTION_PROBLEM_CXX
#endif
