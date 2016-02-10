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
	Nx=20;
	Ny=20;
	Nt=5;


	dx=L/Nx;
	dy=H/Ny;
	
	
	Velocity v(Nx,Ny,L,H);
	//v.ChampsCirculaire(L/2.0,H/2, 10.0);
	//v.ChampsUniformeVx(1.0);
	v.ChampsUniforme(0.5);
	v.WriteGnuPlot("velocity.dat");
	// plot "velocity.dat" u 1:2:3:4 w vec
	cout << "velocity initialise!" << endl;
	//CFL
	dt=(dx/v.max());
	cout << "Vmax = " << v.max() << endl;
	



	Particle n(Nx,Ny,L,H);
	n.InitialSquare(L/2.0,H/2.0,1.0);
	cout << "Particule initialise" << endl;
	n.WriteGnuPlot("particleinit.dat");

	UpwindDCtest1 test1(Nx,Ny,Nt,L,H,tfinal,v,n);


	//LowStorageRungeKuttaIterator timescheme;
	ExplicitEulerIterator timescheme;
	timescheme.SetInitialCondition(0,dt,test1.GetP().Getn(),test1);





	for(int i=0; i<Nt ; i++)
	{
		tn=i*dt;
		timescheme.Advance(i, tn, test1);

	}
	
	n.Setn(timescheme.GetIterate());
	n.WriteGnuPlot("particlefinal.dat");

	return 0;
}
