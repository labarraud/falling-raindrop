#ifndef SPACESCHEME_HXX


#include "define.hxx"
#include "Matrix.hxx"
#include "velocity.hxx"
#include "DiffusionConvectionProblem.hxx"


class UpwindDCtest1 : public DiffusionConvectionProblem

{
public : 
	UpwindDCtest1();
	UpwindDCtest1(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Density& n);
	precision UpwindY(precision dt, precision a, int i, int j, const Matrix& u);
	precision SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u);
  void AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y, const Matrix& sec_membre);
	precision computedt(precision cfl);
};

class UpwindDCOrder2 : public DiffusionConvectionProblem
{
public :
	UpwindDCOrder2();
	UpwindDCOrder2(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Density& n);
	precision UpwindY(precision dt, precision b, int i, int j, const Matrix& u);
	precision SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u);
	void AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y, const Matrix& sec_membre);
	precision computedt(precision cfl);
};

class UpwindDCOrder3 : public DiffusionConvectionProblem

{
public :
	UpwindDCOrder3();
	UpwindDCOrder3(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Density& n);
		precision UpwindY(precision dt, precision a, int i, int j, const Matrix& u);
		precision SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u);
	  void AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y, const Matrix& sec_membre);
	 	precision computedt(precision cfl);
};

class UpwindDCOrder4 : public DiffusionConvectionProblem

{
protected:
public :
	UpwindDCOrder4();
	UpwindDCOrder4(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,precision _D,Density& n);
		precision UpwindY(precision dt, precision a, int i, int j, const Matrix& u);
		precision SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u);
	  void AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y, const Matrix& sec_membre);
	 	precision computedt(precision cfl);
};

class LaxWendroff : public DiffusionConvectionProblem

{
public :
	LaxWendroff();
	LaxWendroff(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Density& n);
	  void AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y, const Matrix& sec_membre);
	 	precision computedt(precision cfl);
};


#define SPACESCHEME_HXX
#endif
