#ifndef SPACESCHEME_CXX

#include "spacescheme.hxx"

UpwindDCtest1::UpwindDCtest1(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Particle& n)
	:	DiffusionConvectionProblem(Nx,Ny,Nt,L,H,tfinal,V,n)
{ }

precision UpwindDCtest1::UpwindY(precision dt, precision b, int i, int j, const Matrix& u)
{
	precision theta(dt/Delta_y);
	int sign_b((b < 0) ? -1 : 1);
	return (u(i,j) - b*theta*sign_b*(u(i,j)-u((Ny+i-sign_b)%Ny,j)));
}

precision UpwindDCtest1 :: SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u)
{
	double sigma(dt/Delta_x), uij;
	int sign_a((a < 0) ? -1 : 1);
	uij = UpwindY(dt, b, i, j, u);
	return (uij - a*sigma*sign_a*(uij-UpwindY(dt, b, i, (Nx+j-sign_a)%Nx, u)));
}

void UpwindDCtest1 :: AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y)
{

	//alpha = delta_t
	precision /*v[4], sigma(alpha/Delta_x), theta(alpha/Delta_y), */a, b;
    for (int i=0; i<Ny; i++)
    {
      for (int j=0; j<Nx; j++)
	{
	  a = velocity.GetVX(i,j);
	  b = velocity.GetVY(i,j);

	  y(i,j) += SplittingX(alpha, a, b, i, j, u) - u(i,j);
/*
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
				+ alpha*D*((u((Ny+i+1)%Ny)(j)-2*u(i)(j)-u((Ny+i-1)%Ny)(j))/(Delta_x*Delta_x)
					+ (u(i)((Nx+j+1)%Nx)-2*u(i)(j)-u(i)((Nx+j-1)%Nx))/(Delta_y*Delta_y));
*/
	}
    }

}

//---------------------------------------------------------------------
// Order 2
//---------------------------------------------------------------------
UpwindDCOrder2::UpwindDCOrder2(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Particle& n)
	:	DiffusionConvectionProblem(Nx,Ny,Nt,L,H,tfinal,V,n)
{ }

precision UpwindDCOrder2::UpwindY(precision dt, precision b, int i, int j, const Matrix& u)
{
	precision theta(dt/Delta_y);
	int sign_b((b < 0) ? -1 : 1);
	return (u(i,j) - b*theta*sign_b*(3.0*u(i,j)-4.0*u((Ny+i-sign_b)%Ny,j)+u((Ny+i-2*sign_b)%Ny,j))/2.0);
}

precision UpwindDCOrder2::SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u)
{
	double sigma(dt/Delta_x), uij;
	int sign_a((a < 0) ? -1 : 1);
	uij = UpwindY(dt, b, i, j, u);
	return (uij - a*sigma*sign_a*(3.0*uij-4.0*UpwindY(dt, b, i, (Nx+j-sign_a)%Nx, u)+UpwindY(dt, b, i, (Nx+j-2*sign_a)%Nx, u))/2.0);
}

void UpwindDCOrder2::AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y)
{

	//alpha = delta_t
    for (int i=0; i<Ny; i++) {
    	for (int j=0; j<Nx; j++) {
			y(i,j) += SplittingX(alpha, velocity.GetVX(i,j), velocity.GetVY(i,j), i, j, u) - u(i,j);
		}
    }
}

//---------------------------------------------------------------------
// Order 3
//---------------------------------------------------------------------

UpwindDCOrder3::UpwindDCOrder3(int Nx,int Ny,int Nt,double L,double H,double tfinal,Velocity& V,Particle& n)
	:	DiffusionConvectionProblem(Nx,Ny,Nt,L,H,tfinal,V,n)
{ }

precision UpwindDCOrder3::UpwindY(double dt, double b, int i, int j, const Matrix& u)
{
	precision theta(dt/Delta_y);
	int sign_b((b < 0) ? -1 : 1);
	return (u(i,j) - b*theta*sign_b*(2.0*u((Ny+i+sign_b)%Ny,j)+3.0*u(i,j)-6.0*u((Ny+i-sign_b)%Ny,j)+u((Ny+i-2*sign_b)%Ny,j))/6.0);
}

precision UpwindDCOrder3::SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u)
{
	double sigma(dt/Delta_x), uij;
	int sign_a((a < 0) ? -1 : 1);
	uij = UpwindY(dt, b, i, j, u);
	return (uij - a*sigma*sign_a*(2.0*UpwindY(dt, b, i, (Nx+j+sign_a)%Nx, u)+3.0*uij-6.0*UpwindY(dt, b, i, (Nx+j-sign_a)%Nx, u)+UpwindY(dt, b, i, (Nx+j-2*sign_a)%Nx, u))/6.0);
}

void UpwindDCOrder3::AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y)
{

	//alpha = delta_t
    for (int i=0; i<Ny; i++) {
    	for (int j=0; j<Nx; j++) {
			y(i,j) += SplittingX(alpha, velocity.GetVX(i,j), velocity.GetVY(i,j), i, j, u) - u(i,j);
		}
    }
}

#define SPACESCHEME_CXX
#endif
