#ifndef SPACESCHEME_CXX
#include "spacescheme.hxx"


void UpwindDCtest1 :: AddFunction(double alpha, const Vector<Vector<double> >& u, double t, Vector<Vector<double> >& y)
{
// LIMIT CONDITION NO DEFINE !!

	//alpha = delta_t
	double v[4], sigma(alpha/Delta_x), theta(alpha/Delta_y);
    for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
	{
	  if (VX(i)(j)>0 && VY(i)(j)>0)
	    {
	      v[0] = u(i)(j);
	      v[1] = u(i)((j-1)%Ny);
	      v[2] = u((i-1%)Nx)(j);
	      v[3] = u((i-1)%Nx)((j-1)%Ny);
	    }
	  
	  
	  else (VX(i)(j)<0 && VY(i)(j)>0)
		{
	      v[0] = u(i)(j);
	      v[1] = u(i)((j+1)%Ny);
	      v[2] = u((i-1)%Nx)(j);
	      v[3] = u((i-1)%Nx)((j+1)%Ny);
		}
	
	  
	  else (VX(i)(j)>0 && VY(i)(j)<0)
		{
	      v[0] = u(i)(j);
	      v[1] = u((i-1)%Nx)(j);
	      v[2] = u(i)((j+1)%Ny);
	      v[3] = u((i+1)%Nx)((j-1)%Ny);
		}

	  else (VX(i)(j)<0 && VY(i)(j)<0)
		{
	      v[0] = u(i)(j);
	      v[1] = u(i)((j+1)%Ny);
	      v[2] = u((i+1)%Nx)(j);
	      v[3] = u((i+1)%Nx)((j+1)%Ny);
	    }
	  y(i)(j) += (1.-VX(i)(j)*sigma)*(1.-VY(i)(j)*theta)*v[0] + (1.-VX(i)(j)*sigma)*VY(i)(j)*theta*v[1]
			+ VX(i)(j)*sigma*(1.-VY(i)(j)*theta)*v[2] + VX(i)(j)*sigma*VY(i)(j)*theta*v[3]
			+ alpha*D*((u(i+1)(j)-2*u(i)(j)-u(i-1)(j))/(Delta_x^2) + (u(i)(j+1)-2*u(i)(j)-u(i)(j-1))/(Delta_y^2));
	}
    }

}

#define SPACESCHEME_CXX
#endif
