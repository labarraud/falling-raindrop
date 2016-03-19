#ifndef SPACESCHEME_HXX


#include "define.hxx"
#include "Matrix.hxx"
#include "velocity.hxx"
#include "DiffusionConvectionProblem.hxx"


class UpwindDCtest1 : public DiffusionConvectionProblem

{
public : 

	UpwindDCtest1(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Particle& n);
	precision UpwindY(precision dt, precision a, int i, int j, const Matrix& u);
	precision SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u);
  void AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y);
};

class UpwindDCOrder2 : public DiffusionConvectionProblem
{
public :

	UpwindDCOrder2(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Particle& n);
	precision UpwindY(precision dt, precision b, int i, int j, const Matrix& u);
	precision SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u);
	void AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y);};

class UpwindDCOrder3 : public DiffusionConvectionProblem

{
public :

	UpwindDCOrder3(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Particle& n);
		precision UpwindY(precision dt, precision a, int i, int j, const Matrix& u);
		precision SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u);
	  void AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y);
};


#define SPACESCHEME_HXX
#endif
