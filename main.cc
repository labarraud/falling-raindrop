#include "linalg/Linalg.hxx"

using namespace linalg;

#include "velocity.cxx"
#include "particle.cxx"
#include "DiffusionConvectionProblem.cxx"
#include "spacescheme.cxx"
#include "TimeScheme.cxx"




int main()
{
	double dx,dy,dt,L,H,tn,tfinal;
	int Nx,Ny,Nt;

	L=10;
	H=10;
	Nx=10;
	Ny=10;
	Nt=10;


	dx=L/Nx;
	dy=H/Ny;
	
	
	Velocity v(Nx,Ny,L,H);
	//v.ChampsCirculaire(L/2.0,H/2, 0.5);
	//v.ChampsUniformeVx(0.5);
	v.ChampsUniforme(0.5);
	v.WriteGnuPlot("velocity.dat");
	// plot "velocity.dat" u 1:2:3:4 w vec
	cout << "velocity initialise!";
	//CFL
	dt=0.9*(dx/v.max());
	



	Particle n(Nx,Ny,L,H);
	n.InitialSquare(2.0,2.0,0.);
	cout << "Particule initialise";
	n.WriteGnuPlot("particleinit.dat");

	UpwindDCtest1 test1(Nx,Ny,Nt,L,H,tfinal,v,n);


	LowStorageRungeKuttaIterator timescheme;
	timescheme.SetInitialCondition(0,dt,test1.GetP().Getn(),test1);





	for(int i=0; i<Nt ; i++)
	{
		tn=i*dt;
		timescheme.Advance(0, tn, test1);

	}
	
	n.Setn(timescheme.GetIterate());
	n.WriteGnuPlot("particlefinal.dat");

	return 0;
}
