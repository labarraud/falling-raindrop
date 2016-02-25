#ifndef SPACESCHEME_CXX
#include "spacescheme.hxx"

UpwindDCtest1::UpwindDCtest1(int Nx,int Ny,int Nt,double L,double H,double tfinal,Velocity& V,Particle& n)
	:	DiffusionConvectionProblem(Nx,Ny,Nt,L,H,tfinal,V,n)
{ }

void UpwindDCtest1 :: AddFunction(double alpha, const Vector<Vector<double> >& u, double t, Vector<Vector<double> >& y)
{

	//alpha = delta_t
	double v[4], sigma(alpha/Delta_x), theta(alpha/Delta_y), a, b;
    for (int i=0; i<Ny; i++)
    {
      for (int j=0; j<Nx; j++)
	{
	  a = velocity.GetVX(i,j);
	  b = velocity.GetVY(i,j);

	  if (a>=0 && b>=0)
	    {
	      v[0] = u(i)(j);
	      v[1] = u((Ny+i-1)%Ny)(j);
	      v[2] = u(i)((Nx+j-1)%Nx);
	      v[3] = u((Ny+i-1)%Ny)((Nx+j-1)%Nx);
	    }
	  
	  
	  else if(a<=0 && b>=0)
		{
	      v[0] = u(i)(j);
	      v[1] = u((Ny+i-1)%Ny)(j);
	      v[2] = u(i)((Nx+j+1)%Nx);
	      v[3] = u((Ny+i-1)%Ny)((Nx+j+1)%Nx);
	      a = -a;
		}
	
	  
	  else if(a>=0 && b<=0)
		{
	      v[0] = u(i)(j);
	      v[1] = u((Ny+i+1)%Ny)(j);
	      v[2] = u(i)((Nx+j-1)%Nx);
	      v[3] = u((Ny+i+1)%Ny)((Nx+j-1)%Nx);
	      b = -b;
		}

	  else// if(a<=0 && b<=0)
		{
	      v[0] = u(i)(j);
	      v[1] = u((Ny+i+1)%Ny)(j);
	      v[2] = u(i)((Nx+j+1)%Nx);
	      v[3] = u((Ny+i+1)%Ny)((Nx+j+1)%Nx);
	      a = -a;
	      b = -b;
	    }
	  y(i)(j) += (1.-a*sigma)*(1.-b*theta)*v[0]-v[0] + (1.-a*sigma)*b*theta*v[1]
			+ a*sigma*(1.-b*theta)*v[2] + a*sigma*b*theta*v[3]
				/*+ alpha*D*((u((Ny+i+1)%Ny)(j)-2*u(i)(j)-u((Ny+i-1)%Ny)(j))/(Delta_x*Delta_x)
					+ (u(i)((Nx+j+1)%Nx)-2*u(i)(j)-u(i)((Nx+j-1)%Nx))/(Delta_y*Delta_y))*/;
	}
    }

}

#define SPACESCHEME_CXX
#endif
