#include "linalg/Linalg.hxx"


#include "velocity.cxx"
#include "particle.cxx"

using namespace linalg;


int main()
{
	double dx,dy,dt,L,H,tfinal;
	int Nx,Ny,Nt;

	L=10;
	H=10;
	Nx=10;
	Ny=10;

	dx=L/Nx;
	dy=H/Ny;
	dt=tfinal/Nt;

	Velocity v(Nx,Ny,L,H);
	v.ChampsCirculaire(L/2.0,H/2, 0.5);
	v.WriteGnuPlot("velocity.dat");
	// plot "velocity.dat" u 1:2:3:4 w vec
	cout << "ok cela tourne!";

	Particle n(Nx,Ny,L,H);
	cout << "ok cela tourne!";
	n.InitialSquare(2.0,2.0,0.5);
	cout << "ok cela tourne!";
	n.WriteGnuPlot("particleinit.dat");





	for(int i=0; i<Nt ; i++)
	{


	}

	return 0;
}
