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

	L=20;
	H=20;
	Nx=200;
	Ny=200;
	Nt=1000;


	dx=L/Nx;
	dy=H/Ny;
	
	
	Velocity v(Nx,Ny,L,H);
	v.ChampsCirculaire(L/2.0,H/2, 1.0);
	//v.ChampsUniformeVx(-1.0);
	//v.ChampsUniforme(-0.5);
	v.WriteGnuPlot("velocity.dat");
	// plot "velocity.dat" u 1:2:3:4 w vec
	cout << "velocity initialise!" << endl;
	//CFL
	dt=(dx/v.max());
	cout << "Vmax = " << v.max() << endl;
	



	Particle n(Nx,Ny,L,H);
	n.InitialSquare(L/3.0,H/3.0,0.5);
	cout << "Particule initialise" << endl;
	n.WriteGnuPlot("particleinit.dat");

	UpwindDCtest1 test1(Nx,Ny,Nt,L,H,tfinal,v,n);


	//LowStorageRungeKuttaIterator timescheme;
	ExplicitEulerIterator timescheme;
	timescheme.SetInitialCondition(0,dt,test1.GetP().Getn(),test1);




	int nDisplay(10);
	for(int i=0; i<Nt ; i++)
	{
		tn=i*dt;
		timescheme.Advance(i, tn, test1);
		if((i%nDisplay)==0) {
			n.Setn(timescheme.GetIterate());
			n.WriteGnuPlot("animate/particle" + to_string(i/nDisplay) + ".dat");
		}

	}
	
	n.Setn(timescheme.GetIterate());
	n.WriteGnuPlot("particlefinal.dat");

	return 0;
}
