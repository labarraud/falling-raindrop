#ifndef SPACESCHEME_CXX
#include "spacescheme.hxx"


void UpwindDCtest1 :: AddFunction(double alpha, const Vector<Vector<double> >& u, double t, Vector<Vector<double> >& y)
{
// LIMIT CONDITION NO DEFINE !!

	//alpha = delta_t
	double v[4], sigma(alpha/Delta_x), theta(alpha/Delta_y), a, b;
    for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
	{
	  a = velocity.GetVX(i,j);
	  b = velocity.GetVY(i,j);

	  if (a>0 && b>0)
	    {
	      v[0] = u(i)(j);
	      v[1] = u(i)((j-1)%Ny);
	      v[2] = u((i-1%)Nx)(j);
	      v[3] = u((i-1)%Nx)((j-1)%Ny);
	    }
	  
	  
	  else (a<0 && b>0)
		{
	      v[0] = u(i)(j);
	      v[1] = u(i)((j+1)%Ny);
	      v[2] = u((i-1)%Nx)(j);
	      v[3] = u((i-1)%Nx)((j+1)%Ny);
	      a = -a;
		}
	
	  
	  else (a>0 && b<0)
		{
	      v[0] = u(i)(j);
	      v[1] = u((i-1)%Nx)(j);
	      v[2] = u(i)((j+1)%Ny);
	      v[3] = u((i+1)%Nx)((j-1)%Ny);
	      b = -b;
		}

	  else (a<0 && b<0)
		{
	      v[0] = u(i)(j);
	      v[1] = u(i)((j+1)%Ny);
	      v[2] = u((i+1)%Nx)(j);
	      v[3] = u((i+1)%Nx)((j+1)%Ny);
	      a = -a;
	      b = -b;
	    }
	  y(i)(j) += (1.-a*sigma)*(1.-b*theta)*v[0] + (1.-a*sigma)*b*theta*v[1]
			+ a*sigma*(1.-b*theta)*v[2] + a*sigma*b*theta*v[3]
			+ alpha*D*((u(i+1)(j)-2*u(i)(j)-u(i-1)(j))/(Delta_x^2) + (u(i)(j+1)-2*u(i)(j)-u(i)(j-1))/(Delta_y^2));
	}
    }

}

#define SPACESCHEME_CXX
#endif
