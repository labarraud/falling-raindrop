#ifndef SPACESCHEME_HXX



class DecentrerDC : public DiffusionConvectionProblem
{
public : 

  virtual void AddFunction(double alpha, const Vector<Vector<double>>& u, double t, Vector<Vector<double>>& y);
}


#define SPACESCHEME_HXX
#endif
