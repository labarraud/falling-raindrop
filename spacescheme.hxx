#ifndef SPACESCHEME_HXX



class UpwindDCtest1 : public DiffusionConvectionProblem

{
public : 

	UpwindDCtest1(int Nx,int Ny,int Nt,double L,double H,double tfinal,Velocity& V,Particle& n);
  void AddFunction(double alpha, const Vector<Vector<double> >& u, double t, Vector<Vector<double> >& y);
};




#define SPACESCHEME_HXX
#endif
