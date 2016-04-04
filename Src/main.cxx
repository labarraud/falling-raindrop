#include <string>

#include "../Include/Matrix.hxx"
#include "../Include/velocity.hxx"
#include "../Include/particle.hxx"
#include "../Include/DiffusionConvectionProblem.hxx"
#include "../Include/spacescheme.hxx"
#include "../Include/TimeScheme.hxx"




int main()
{

/*
 * coucou
  //-----------experimental order -------------

	precision mindxy,hdxy,maxdxy,cfl,tmaxdemi,omega;
	UpwindDCtest1 ode;
	LowStorageRungeKuttaIterator time;

	mindxy=0.0015;
	hdxy=0.005;
	maxdxy=0.09;
	cfl=1;
	tmaxdemi=1;
	omega=5.0;

	error_orderxy_circle(mindxy,hdxy,maxdxy,cfl , tmaxdemi,omega, ode, time, "error_upwind1.dat");
*/
	precision dx,dy,dt,L,H,tn,tfinal,cfl, D;
		int Nx,Ny,Nt;

		L=5;
		H=5;
		Nx=200;
		Ny=200;
		Nt=20000;
		cfl=0.4;

		dx=L/Nx;
		dy=H/Ny;
		D = 0.0018;


		Velocity v(Nx,Ny,L,H);
		v.ChampsCirculaire(L/2.0,H/2, 5.0);
		//v.ChampsUniformeVx(50.0);
		//v.ChampsUniforme(-0.5);
		v.WriteGnuPlot("velocity.dat");
		// plot "velocity.dat" u 1:2:3:4 w vec
		cout << "velocity initialise!" << endl;
		//CFL
		dt=((max(dx,dy)*max(dx,dy)*cfl)/v.max());
		//dt=((max(dx,dy)*max(dx,dy)*cfl)/D);
		tfinal=Nt*dt;
		cout << "Vmax = " << v.max() << endl;



		Density n(Nx,Ny,L,H);
		//n.InitialSquare(L/3.0,H/3.0,0.5);
		n.InitialCircle(L/3.0,H/3.0,0.25);
		//n.InitialGauss(L/3.0,H/3.0,0.5);
		cout << "Particule initialise" << endl;
		n.WriteGnuPlot("particleinit.dat");

		//UpwindDCtest1 test1(Nx,Ny,Nt,L,H,tfinal,v,n);
		//UpwindDCOrder2 test1(Nx,Ny,Nt,L,H,tfinal,v,n);
		//UpwindDCOrder3 test1(Nx,Ny,Nt,L,H,tfinal,v,n);
		UpwindDCOrder4 test1(Nx,Ny,Nt,L,H,tfinal,v,D,n);
		//LaxWendroff test1(Nx,Ny,Nt,L,H,tfinal,v,n);


		LowStorageRungeKuttaIterator timescheme;
		//ExplicitEulerIterator timescheme;
		timescheme.SetInitialCondition(0,dt,test1.GetP(),test1);




		int nDisplay(100);
		/*string file = "scriptan.gnuplot";
		ofstream file_out(file.data());
		file_out.precision(15);*/
		string var,var2;
		for(int i=0; i<Nt ; i++)
		{
			tn=i*dt;
			timescheme.Advance(i, tn, test1);
			if((i%nDisplay)==0) {
				var=(i/nDisplay < 10 ? "0" : "");
				var2=(i/nDisplay < 100 ? "0" : "");
				/*file_out << "set terminal postscript eps enhanced color" << endl;
				file_out << "set output '" << ("animate/particle" + var + var2 + to_string(i/nDisplay) + ".eps'") << endl;
				file_out << "set pm3d map" << endl;
				file_out << ("splot 'animate/particle" + to_string(i/nDisplay) + ".dat' matrix") << endl  << endl;*/

				//static_cast<const Density&>(timescheme.GetIterate()).WriteGnuPlot("animate/particle" + to_string(i/nDisplay) + ".dat");

				static_cast<const Density&>(timescheme.GetIterate()).WriteVtk("vtk/particle" + var + var2 + to_string(i/nDisplay) + ".vtk");
			}

		}

		/*
		Velocity v2(Nx,Ny,L,H);
		v2.ChampsCirculaire(L/2.0,H/2, -5.0);
		//v.ChampsUniformeVx(-1.0);
		//v.ChampsUniforme(-0.5);
		v2.WriteGnuPlot("velocity.dat");
		// plot "velocity.dat" u 1:2:3:4 w vec
		cout << "velocity initialise!" << endl;
		//CFL
		dt=(dx/v.max());
		cout << "Vmax = " << v2.max() << endl;

		UpwindDCtest1 test2(Nx,Ny,Nt,L,H,tfinal,v2,n);
		//LowStorageRungeKuttaIterator timescheme;
		ExplicitEulerIterator timescheme2;
		timescheme2.SetInitialCondition(0,dt,test2.GetP().Getn(),test2);
		for(int i=Nt; i<2*Nt ; i++)
		{
			tn=i*dt;
			timescheme2.Advance(i, tn, test2);
			if((i%nDisplay)==0) {
				var=(i/nDisplay < 10 ? "0" : "");
				var2=(i/nDisplay < 100 ? "0" : "");
				file_out << "set terminal postscript eps enhanced color" << endl;
				file_out << "set output '" << ("animate/particle" + var + var + to_string(i/nDisplay) + ".eps'") << endl;
				file_out << "set pm3d map" << endl;
				file_out << ("splot 'animate/particle" + to_string(i/nDisplay) + ".dat' matrix") << endl  << endl;
				n.Setn(timescheme2.GetIterate());
				n.WriteGnuPlot("animate/particle" + to_string(i/nDisplay) + ".dat");
			}
		}
		 file_out.close();*/
		//n.Setn(timescheme.GetIterate());
		n.WriteGnuPlot("particlefinal.dat");

		return 0;
}
