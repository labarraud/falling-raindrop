#ifndef SPACESCHEME_CXX
#include "spacescheme.hxx"


void DecentrerDC :: AddFunction(double alpha, const Vector<Vector<double>>& u, double t, Vector<Vector<double>>& y)
{

    for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
	{
	  if (VX(i)(j)>0 && VY(i)(j)>0)
	    {
	      y(i,j) +=alpha*(VX(i)(j)*(u(i)(j)-u(i-1)(j))/Delta_x +  VY(i)(j)*(u(i)(j)-u(i)(j-1))/Delta_y - D*( (u(i+1)(j)-2*u(i)(j)-u(i-1)(j))/((Delta_x)^2) +  (u(i)(j+1)-2*u(i)(j)-u(i)(j-1))/((Delta_y)^2)));
	    }
	  
	  
	  else (VX(i)(j)<0 && VY(i)(j)>0)
		{
		  y(i,j) += alpha*(VX(i)(j)*(u(i+1)(j)-u(i)(j))/Delta_x +  VY(i)(j)*(u(i)(j)-u(i)(j-1))/Delta_y - D*( (u(i+1)(j)-2*u(i)(j)-u(i-1)(j))/((Delta_x)^2) +  (u(i)(j+1)-2*u(i)(j)-u(i)(j-1))/((Delta_y)^2)));
		}
	
	  
	  else (VX(i)(j)>0 && VY(i)(j)<0)
		{
		  y(i,j) += alpha*(VX(i)(j)*(u(i)(j)-u(i-1)(j))/Delta_x +  VY(i)(j)*(u(i)(j+1)-u(i)(j))/Delta_y - D*( (u(i+1)(j)-2*u(i)(j)-u(i-1)(j))/((Delta_x)^2) +  (u(i)(j+1)-2*u(i)(j)-u(i)(j-1))/((Delta_y)^2)));
		}

	  else (VX(i)(j)<0 && VY(i)(j)<0)
		{
		  y(i,j) += alpha*(VX(i)(j)*(u(i+1)(j)-u(i)(j))/Delta_x +  VY(i)(j)*(u(i)(j+1)-u(i)(j))/Delta_y - D*( (u(i+1)(j)-2*u(i)(j)-u(i-1)(j))/((Delta_x)^2) +  (u(i)(j+1)-2*u(i)(j)-u(i)(j-1))/((Delta_y)^2)));
	    }
	}
    }

}

#define SPACESCHEME_CXX
#endif
