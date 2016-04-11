#ifndef SPACESCHEME_CXX

#include "../Include/spacescheme.hxx"

UpwindDCtest1::UpwindDCtest1() {}

UpwindDCtest1::UpwindDCtest1(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Density& n)
	:	DiffusionConvectionProblem(Nx,Ny,Nt,L,H,tfinal,V,n)
{ }

precision UpwindDCtest1::UpwindY(precision dt, precision b, int i, int j, const Matrix& u)
{
	precision theta(dt/dy);
	int sign_b((b < 0) ? -1 : 1);
	return (u(i,j) - b*theta*sign_b*(u(i,j)-u((Ny+i-sign_b)%Ny,j)));
}

precision UpwindDCtest1 :: SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u)
{
	double sigma(dt/dx), uij;
	int sign_a((a < 0) ? -1 : 1);
	uij = UpwindY(dt, b, i, j, u);
	return (uij - a*sigma*sign_a*(uij-UpwindY(dt, b, i, (Nx+j-sign_a)%Nx, u)));
}

void UpwindDCtest1 :: AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y)
{

	//alpha = delta_t
	precision v[4], sigma(alpha/dx), theta(alpha/dy), a, b;
    for (int i=0; i<Ny; i++)
    {
      for (int j=0; j<Nx; j++)
	{
	  a = velocity.GetVX(i,j);
	  b = velocity.GetVY(i,j);

	  //y(i,j) += SplittingX(alpha, a, b, i, j, u) - u(i,j);

	  if (a>=0 && b>=0)
	    {
	      v[0] = u(i,j);
	      v[1] = u((Ny+i-1)%Ny,j);
	      v[2] = u(i,(Nx+j-1)%Nx);
	      v[3] = u((Ny+i-1)%Ny,(Nx+j-1)%Nx);
	    }
	  
	  
	  else if(a<=0 && b>=0)
		{
	      v[0] = u(i,j);
	      v[1] = u((Ny+i-1)%Ny,j);
	      v[2] = u(i,(Nx+j+1)%Nx);
	      v[3] = u((Ny+i-1)%Ny,(Nx+j+1)%Nx);
	      a = -a;
		}
	
	  
	  else if(a>=0 && b<=0)
		{
	      v[0] = u(i,j);
	      v[1] = u((Ny+i+1)%Ny,j);
	      v[2] = u(i,(Nx+j-1)%Nx);
	      v[3] = u((Ny+i+1)%Ny,(Nx+j-1)%Nx);
	      b = -b;
		}

	  else// if(a<=0 && b<=0)
		{
	      v[0] = u(i,j);
	      v[1] = u((Ny+i+1)%Ny,j);
	      v[2] = u(i,(Nx+j+1)%Nx);
	      v[3] = u((Ny+i+1)%Ny,(Nx+j+1)%Nx);
	      a = -a;
	      b = -b;
	    }
	  y(i,j) += (1.-a*sigma)*(1.-b*theta)*v[0]-v[0] + (1.-a*sigma)*b*theta*v[1]
			+ a*sigma*(1.-b*theta)*v[2] + a*sigma*b*theta*v[3]
				+ alpha*D*((u((Ny+i+1)%Ny,j)-2*u(i,j)-u((Ny+i-1)%Ny,j))/(dx*dx)
					+ (u(i,(Nx+j+1)%Nx)-2*u(i,j)-u(i,(Nx+j-1)%Nx))/(dy*dy));

	}
    }

}




precision UpwindDCtest1::computedt(precision cfl){

	return ((max(dx,dy)*cfl)/velocity.max());
}


//---------------------------------------------------------------------
// Order 2
//---------------------------------------------------------------------
UpwindDCOrder2::UpwindDCOrder2() {}

UpwindDCOrder2::UpwindDCOrder2(int Nx,int Ny,int Nt,precision L,precision H,precision tfinal,Velocity& V,Density& n)
	:	DiffusionConvectionProblem(Nx,Ny,Nt,L,H,tfinal,V,n)
{ }

precision UpwindDCOrder2::UpwindY(precision dt, precision b, int i, int j, const Matrix& u)
{
	precision theta(dt/dy);
	int sign_b((b < 0) ? -1 : 1);
	return (u(i,j) - b*theta*sign_b*(3.0*u(i,j)-4.0*u((Ny+i-sign_b)%Ny,j)+u((Ny+i-2*sign_b)%Ny,j))/2.0);
}

precision UpwindDCOrder2::SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u)
{
	double sigma(dt/dx), uij;
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


precision UpwindDCOrder2::computedt(precision cfl){

	return ((max(dx,dy)*max(dx,dy)*cfl)/velocity.max());
}


//---------------------------------------------------------------------
// Order 3
//---------------------------------------------------------------------
UpwindDCOrder3::UpwindDCOrder3() { }
UpwindDCOrder3::UpwindDCOrder3(int Nx,int Ny,int Nt,double L,double H,double tfinal,Velocity& V,Density& n)
	:	DiffusionConvectionProblem(Nx,Ny,Nt,L,H,tfinal,V,n)
{ }

precision UpwindDCOrder3::UpwindY(double dt, double b, int i, int j, const Matrix& u)
{
	precision theta(dt/dy);
	int sign_b((b < 0) ? -1 : 1);
	return (u(i,j) - b*theta*sign_b*(2.0*u((Ny+i+sign_b)%Ny,j)+3.0*u(i,j)-6.0*u((Ny+i-sign_b)%Ny,j)+u((Ny+i-2*sign_b)%Ny,j))/6.0);
}

precision UpwindDCOrder3::SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u)
{
	double sigma(dt/dx), uij;
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

precision UpwindDCOrder3::computedt(precision cfl){

	return ((max(dx,dy)*max(dx,dy)*max(dx,dy)*cfl)/velocity.max());
}

//---------------------------------------------------------------------
// Order 4
//---------------------------------------------------------------------
UpwindDCOrder4::UpwindDCOrder4() { }
UpwindDCOrder4::UpwindDCOrder4(int Nx,int Ny,int Nt,double L,double H,double tfinal,Velocity& V,precision _D,Density& n)
	:	DiffusionConvectionProblem(Nx,Ny,Nt,L,H,tfinal,V,n)
{
	D = _D;
}

precision UpwindDCOrder4::UpwindY(double dt, double b, int i, int j, const Matrix& u)
{
	precision theta(dt/dy);
	int sign_b((b < 0) ? -1 : 1);
	return (u(i,j) - b*theta*sign_b*(3.0*u((Ny+i+sign_b)%Ny,j)+10.0*u(i,j)-18.0*u((Ny+i-sign_b)%Ny,j)+6.0*u((Ny+i-2*sign_b)%Ny,j)-u((Ny+i-3*sign_b)%Ny,j))/12.0);
}

precision UpwindDCOrder4::SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u)
{
	double sigma(dt/dx), uij;
	int sign_a((a < 0) ? -1 : 1);
	uij = UpwindY(dt, b, i, j, u);
	return (uij - a*sigma*sign_a*(3.0*UpwindY(dt, b, i, (Nx+j+sign_a)%Nx, u)+10.0*uij-18.0*UpwindY(dt, b, i, (Nx+j-sign_a)%Nx, u)+6.0*UpwindY(dt, b, i, (Nx+j-2*sign_a)%Nx, u)-UpwindY(dt, b, i, (Nx+j-3*sign_a)%Nx, u))/12.0);
}

void UpwindDCOrder4::AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y)
{

	//alpha = delta_t
    for (int i=0; i<Ny; i++) {
    	for (int j=0; j<Nx; j++) {
			y(i,j) += SplittingX(alpha, velocity.GetVX(i,j), velocity.GetVY(i,j), i, j, u) - u(i,j)
						+ alpha*D*(16.0*u((Ny+i+1)%Ny,j)-u((Ny+i+2)%Ny,j)-30.0*u(i,j)-u((Ny+i-2)%Ny,j)+16.0*u((Ny+i-1)%Ny,j))/(dx*dx)
						+ alpha*D*(16.0*u(i,(Nx+j+1)%Nx)-u(i,(Nx+j+2)%Nx)-30.0*u(i,j)-u(i,(Nx+j-2)%Nx)+16.0*u(i,(Nx+j-1)%Nx))/(dy*dy);
		}
   }
}

precision UpwindDCOrder4::computedt(precision cfl){

	return (max(dx,dy)*(max(dx,dy)*max(dx,dy)*max(dx,dy)*cfl)/velocity.max());
}

//---------------------------------------------------------------------
// LaxWendroff
//---------------------------------------------------------------------
LaxWendroff::LaxWendroff(){ }

LaxWendroff::LaxWendroff(int Nx,int Ny,int Nt,double L,double H,double tfinal,Velocity& V,Density& n)
	:	DiffusionConvectionProblem(Nx,Ny,Nt,L,H,tfinal,V,n)
{ }

void LaxWendroff::AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y)
{

	//alpha = delta_t
	precision sigma(alpha/dx), theta(alpha/dy), a, b;
    for (int i=0; i<Ny; i++) {
    	for (int j=0; j<Nx; j++) {
    		a = velocity.GetVX(i,j);
    		b = velocity.GetVY(i,j);
			y(i,j) += -0.5*sigma*b*(u((Ny+i+1)%Ny,j)-u((Ny+i-1)%Ny,j))
					- 0.5*theta*a*(u(i,(Nx+j+1)%Nx)-u(i,(Nx+j-1)%Nx))
					+ 0.5*sigma*sigma*b*b*(u((Ny+i+1)%Ny,j)+u((Ny+i-1)%Ny,j)-2.0*u(i,j))
					+ 0.5*theta*theta*a*a*(u(i,(Nx+j+1)%Nx)+u(i,(Nx+j-1)%Nx)-2.0*u(i,j))
					+ 0.25*sigma*theta*a*b*((u((Ny+i+1)%Ny,(Nx+j+1)%Nx)-u((Ny+i-1)%Ny,(Nx+j+1)%Nx))-(u((Ny+i+1)%Ny,(Nx+j-1)%Nx)-u((Ny+i-1)%Ny,(Nx+j-1)%Nx)));

    	}
    }
	//cout << "norme de y = " << y.norme2() << endl;
}

precision LaxWendroff::computedt(precision cfl){

	return ((max(dx,dy)*max(dx,dy)*cfl)/velocity.max());
}

#define SPACESCHEME_CXX
#endif
