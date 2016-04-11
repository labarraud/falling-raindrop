#ifndef DIFFUSION_CONVECTION_PROBLEM_CXX

#include "../Include/DiffusionConvectionProblem.hxx"


//! Destructeur
VirtualOdeSystem::~VirtualOdeSystem()
{
}




/*	Vector<double> step_x, step_y;
	Velocity V;
	Density density;
	int Nx,Ny,Nt;
	double L,H,D,dx,dy,dt;*/

DiffusionConvectionProblem::DiffusionConvectionProblem(){ }

	DiffusionConvectionProblem::DiffusionConvectionProblem(int Nx,int Ny,int Nt,double L,double H,double tfinal,Velocity& V,Density& n)
		:	velocity(V), density(n)
	{

		this->Nx=Nx;
		this->Ny=Ny;
		this->L=L;
    	this->H=H;
    	this->D=0;
    	this->tfinal=H;
	    this->dx=L/Nx;
		this->dy=H/Ny;
		this->dt=tfinal/Nt;


	}

	void DiffusionConvectionProblem::SetInitialCondition(int Nx,int Ny,int Nt,double L,double H,double tfinal,Velocity& V,Density& n)
	{
			this->Nx=Nx;
			this->Ny=Ny;
			this->L=L;
	    	this->H=H;
	    	this->D=0;
	    	this->tfinal=H;
		    this->dx=L/Nx;
			this->dy=H/Ny;
			this->dt=tfinal/Nt;
			this->velocity=V;
			this->density=n;
	}

	DiffusionConvectionProblem::~DiffusionConvectionProblem()
	{

	}

	/*void DiffusionConvectionProblem::Init(int Nx, int Ny,double D)
	{
		step_x.Reallocate(Nx+1);
		step_y.Reallocate(Ny+1);
		dx = D/Nx;
		dy = D/Ny;
		L = H = this.D = D;
		if(V != NULL) {
			delete V;
		}
		V = new Velocity(Nx, Ny, L, H);
		for(int i(0); i < Nx+1; ++i) {
			step_x(i) = i*dx;
		}
		for(int i(0); i < Ny+1; ++i) {
			step_y(i) = i*dy;
		}
	}
*/

	void DiffusionConvectionProblem::WriteGnuPlot(const Matrix& M,const string& nom)
	{
		ofstream file_out(nom.data());
		file_out.precision(15);
		for(int i(0), j; i < M.GetN(); i++) {
			for(j = 0; j < M.GetM(); ++j) {
				file_out << step_x[(unsigned)j] << " " << step_y[(unsigned)i] << " " << M(i,j) << '\n';
			}
			file_out  << '\n';
		}

		file_out.close();
	}

	double DiffusionConvectionProblem::GetX(int i) const
	{
		return step_x[(unsigned)i];
	}

	double DiffusionConvectionProblem::GetY(int i) const
	{
		return step_y[(unsigned)i];
	}
	
	Density & DiffusionConvectionProblem::GetP()
	{
		return density;
		
	}

	Velocity & DiffusionConvectionProblem::GetV()

	{
			return velocity;

	}

#define DIFFUSION_CONVECTION_PROBLEM_CXX
#endif
