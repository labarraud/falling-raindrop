#ifndef SPACESCHEME_HXX



class UpwindDCtest1 : public DiffusionConvectionProblem

{
public : 

	UpwindDCtest1(int Nx,int Ny,int Nt,double L,double H,double tfinal,Velocity& V,Particle& n);
	double UpwindY(double dt, double a, int i, int j, const Vector<Vector<double> >& u);
	double SplittingX(double dt, double a, double b, int i, int j, const Vector<Vector<double> >& u);
  void AddFunction(double alpha, const Vector<Vector<double> >& u, double t, Vector<Vector<double> >& y);
};

class UpwindDCOrder2 : public DiffusionConvectionProblem

{
public :

	UpwindDCOrder2(int Nx,int Ny,int Nt,double L,double H,double tfinal,Velocity& V,Particle& n);
	double UpwindY(double dt, double a, int i, int j, const Vector<Vector<double> >& u);
	double SplittingX(double dt, double a, double b, int i, int j, const Vector<Vector<double> >& u);
  void AddFunction(double alpha, const Vector<Vector<double> >& u, double t, Vector<Vector<double> >& y);
};

class UpwindDCOrder3 : public DiffusionConvectionProblem

{
public :

	UpwindDCOrder3(int Nx,int Ny,int Nt,double L,double H,double tfinal,Velocity& V,Particle& n);
	double UpwindY(double dt, double a, int i, int j, const Vector<Vector<double> >& u);
	double SplittingX(double dt, double a, double b, int i, int j, const Vector<Vector<double> >& u);
  void AddFunction(double alpha, const Vector<Vector<double> >& u, double t, Vector<Vector<double> >& y);
};


#define SPACESCHEME_HXX
#endif
