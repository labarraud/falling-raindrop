#ifndef SPACESCHEME_HXX



class UpwindDCtest1 : public DiffusionConvectionProblem

{
public : 

  virtual void AddFunction(double alpha, const Vector<Vector<double> >& u, double t, Vector<Vector<double> >& y);
}


#define SPACESCHEME_HXX
#endif
