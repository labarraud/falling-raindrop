#ifndef DIFFUSION_CONVECTION_PROBLEM_HXX

#include <vector>

#include "define.hxx"
#include "Matrix.hxx"
#include "particle.hxx"
#include "velocity.hxx"


//! Classe abstraite pour la definition d'une EDO
/*!
  une EDO s'ecrit rho' = f(t, rho)  
 */
class VirtualOdeSystem
{
public:
  virtual ~VirtualOdeSystem();

  virtual void SetInitialCondition(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Density& n)=0;

  // fonction pour calculer y = y + alpha f(t, rho)
  virtual void AddFunction(precision alpha, const Matrix& rho, precision t, Matrix& y, const vector<precision>& sec_membre) = 0;
  
};

class DiffusionConvectionProblem : public VirtualOdeSystem
{

protected:
	vector<precision> step_x, step_y;
	Velocity velocity;
	Density density;
	int Nx,Ny,Nt;
	double L,H,tfinal,D,dx,dy,dt;


public:
	DiffusionConvectionProblem();
	DiffusionConvectionProblem(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Density& n);
	void SetInitialCondition(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Density& n);
	virtual ~DiffusionConvectionProblem();
	//void Init(int Nx, int Ny,double D);
	virtual void AddFunction(precision alpha, const Matrix& rho, precision t, Matrix& y, const vector<precision>& sec_membre) = 0;
	void WriteGnuPlot(const Matrix& M,const string& nom);
	double GetX(int i) const;
	double GetY(int i) const;
	Density& GetP();
	Velocity & GetV();
	virtual precision computedt(precision cfl)=0;
};


#define DIFFUSION_CONVECTION_PROBLEM_HXX
#endif
