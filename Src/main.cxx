#include <string>

#include "../Include/Matrix.hxx"
#include "../Include/velocity.hxx"
#include "../Include/particle.hxx"
#include "../Include/DiffusionConvectionProblem.hxx"
#include "../Include/spacescheme.hxx"
#include "../Include/TimeScheme.hxx"
#include "../Include/NavierStokes.hxx"

void testNS2()
{
		precision dx,dy,dt,L,H,tn,tfinal,cfl;
		int Nx,Ny,Nt;

		L=0.075; // tailles en m
		H=0.15;
		Nx=200;
		Ny=400;
		Nt=50000;
		cfl=0.9;

		dx=L/Nx;
		dy=H/Ny;


		Velocity v(Nx,Ny,L,H);
		v.ChampsCircle(0.025,H-0.003,0.0015,0.0,-5.0);
		
		//CFL
		//dt=((max(dx,dy)*cfl)/v.max());
		dt = 0.00005;
		//cout << "dt = " << dt << endl;
		tfinal=Nt*dt;
		//cout << "Vmax = " << v.max() << endl;



		Density n(Nx,Ny,L,H);		
		precision rho_mer(1032.0), p_atm(1015000.0), y, g(9.81);
		
		Matrix p(Ny+1,Nx+1);
		for(int i(0); i < Ny+1; ++i) {
			y = H-i*dy;
			for(int j(0); j < Nx+1; ++j) {
				n(i,j) = rho_mer;
				p(i,j) = p_atm + rho_mer*g*y;
			}
		}
		NavierStokes2 ns(Nx,Ny,Nt,L,H,tfinal,v,n,p);
		int nDisplay(10);
		string var,var2;
		for(int i=0; i<Nt ; i++)
		{
			cout << "iteration " << i << "  -  p residu = ";
			tn=i*dt;
			//timescheme.Advance(i, tn, test1);
			ns.Advance(i, tn);
			if((i%nDisplay)==0) {
				var=(i/nDisplay < 10 ? "0" : "");
				var2=(i/nDisplay < 100 ? "0" : "");
				ns.WriteVtk("vtk2/particle" + var + var2 + to_string(i/nDisplay) + ".vtk");
			}

		}
}

void testrotateDC()
{

	precision dx,dy,dt,L,H,tn,tfinal,cfl;
	int Nx,Ny,Nt;

	L=10;
	H=10;
	Nx=100;
	Ny=100;
	Nt=20000;
	cfl=	0.4;

	dx=L/Nx;
	dy=H/Ny;


	Velocity v(Nx,Ny,L,H);
	//v.ChampsCirculaire(L/2.0,H/2, 5.0);
	v.ChampsUniformeVx(-1.0);
	//v.ChampsUniforme(-0.5);
	v.WriteGnuPlot("velocity.dat");
	// plot "velocity.dat" u 1:2:3:4 w vec
		cout << "velocity initialise!" << endl;
	//CFL
			dt=((max(dx,dy)*max(dx,dy)*cfl)/v.max());
			tfinal=Nt*dt;
			cout << "Vmax = " << v.max() << endl;



			Density n(Nx,Ny,L,H);
			//n.InitialSquare(L/3.0,H/3.0,0.5);
			n.InitialCircle(L/3.0,H/3.0,0.5);
			//n.InitialGauss(L/3.0,H/3.0,0.5);
			cout << "Particule initialise" << endl;
			n.WriteGnuPlot("particleinit.dat");

			//UpwindDCtest1 test1(Nx,Ny,Nt,L,H,tfinal,v,n);
			//UpwindDCOrder2 test1(Nx,Ny,Nt,L,H,tfinal,v,n);
			//UpwindDCOrder3 test1(Nx,Ny,Nt,L,H,tfinal,v,n);
			LaxWendroff test1(Nx,Ny,Nt,L,H,tfinal,v,n);


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


					//n.WriteGnuPlot("animate/particle" + to_string(i/nDisplay) + ".dat");

					timescheme.GetIterate().WriteVtk("vtk" + var + var2 + to_string(i/nDisplay) + ".vtk", dx, dy);
				}

			}
		}



void testrotatedirDC()
{
	precision dx,dy,dt,L,H,tn,tfinal,cfl;
	int Nx,Ny,Nt;

	L=10;
	H=10;
	Nx=100;
	Ny=100;
	Nt=20000;
	cfl=	0.4;

	dx=L/Nx;
	dy=H/Ny;


	Velocity v(Nx,Ny,L,H);
	//v.ChampsCirculaire(L/2.0,H/2, 5.0);
	v.ChampsUniformeVx(-1.0);
	//v.ChampsUniforme(-0.5);
	v.WriteGnuPlot("velocity.dat");
	// plot "velocity.dat" u 1:2:3:4 w vec
		cout << "velocity initialise!" << endl;
	//CFL
			dt=((max(dx,dy)*max(dx,dy)*cfl)/v.max());
			tfinal=Nt*dt;
			cout << "Vmax = " << v.max() << endl;



			Density n(Nx,Ny,L,H);
			//n.InitialSquare(L/3.0,H/3.0,0.5);
			n.InitialCircle(L/3.0,H/3.0,0.5);
			n.SetBoundaryCondition(dirichlet,0,dirichlet,0,periodic,0,periodic,0);
			//n.InitialGauss(L/3.0,H/3.0,0.5);
			cout << "Particule initialise" << endl;
			n.WriteGnuPlot("particleinit.dat");

			UpwindDCtest1 test1(Nx,Ny,Nt,L,H,tfinal,v,n);
			//UpwindDCOrder2 test1(Nx,Ny,Nt,L,H,tfinal,v,n);
			//UpwindDCOrder3 test1(Nx,Ny,Nt,L,H,tfinal,v,n);
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


					//n.WriteGnuPlot("animate/particle" + to_string(i/nDisplay) + ".dat");

					timescheme.GetIterate().WriteVtk("vtkdirichlet/vtk" + var + var2 + to_string(i/nDisplay) + ".vtk", dx, dy);
				}

			}
		}





int main()
{

/*
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



	testrotateDC();
	//testrotatedirDC();


		return 0;
}
