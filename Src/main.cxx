#include <string>

#include "../Include/Matrix.hxx"
#include "../Include/velocity.hxx"
#include "../Include/particle.hxx"
#include "../Include/DiffusionConvectionProblem.hxx"
#include "../Include/spacescheme.hxx"
#include "../Include/TimeScheme.hxx"
#include "../Include/NavierStokes2.hxx"

void testNS2()
{
		precision dx,dy,dt,L,H,tn,tfinal,cfl;
		int Nx,Ny,Nt;

		L=0.05; // tailles en m
		H=0.1;
		Nx=100;
		Ny=200;
		Nt=10000;
		cfl=0.9;

		dx=L/Nx;
		dy=H/Ny;


		Velocity v(Nx,Ny,L,H);
		precision Xcenter(L/2), Ycenter(H-0.0015), radius(0.0015);
		v.ChampsCircle(Xcenter,Ycenter,radius,0.0,-5.0);
		
		//CFL
		//dt=(max(dx,dy)*(max(dx,dy)*cfl))/v.max();
		dt = 0.00005;
		//cout << "dt = " << dt << endl;
		tfinal=Nt*dt;
		//cout << "Vmax = " << v.max() << endl;



		Density n(Nx,Ny,L,H);
		precision rho_mer(1032.0), p_atm(1015000.0), y, g(9.81);
		
		v.GetAllVX().SetBoundaryCondition(neumann,0,neumann,0,neumann,0,neumann,0);
		v.GetAllVY().SetBoundaryCondition(neumann,0,neumann,0,neumann,0,neumann,0);
		Matrix p(Ny+1,Nx+1);
		p.SetBoundaryCondition(neumann,0,neumann,0,dirichlet,p_atm+rho_mer*g*(H+dy),dirichlet,p_atm-rho_mer*g*dy);
		for(int i(0); i < Ny+1; ++i) {
			y = H-i*dy;
			for(int j(0); j < Nx+1; ++j) {
				n(i,j) = rho_mer;
				p(i,j) = p_atm + rho_mer*g*y;
			}
		}
		NavierStokes2 ns(Nx,Ny,Nt,L,H,tfinal,v,n,p);
		int nDisplay(50);
		string var,var2;
		for(int i=0; i<Nt ; i++)
		{
			cout << "iteration " << i << "  -  p residu = ";
			tn=i*dt;
			ns.Advance(i, tn);
			if((i%nDisplay)==0) {
				var=(i/nDisplay < 10 ? "0" : "");
				var2=(i/nDisplay < 100 ? "0" : "");
				ns.WriteVtk("output_test/testNS2/particle" + var + var2 + to_string(i/nDisplay) + ".vtk");
			}

		}
}

void testpoiseuille()
{
		precision dx,dy,dt,L,H,tn,tfinal,cfl;
		int Nx,Ny,Nt;

		L=0.1; // tailles en m
		H=1;
		Nx=20;
		Ny=100;
		Nt=50000;
		cfl=0.1;

		dx=L/Nx;
		dy=H/Ny;


		Velocity v(Nx,Ny,L,H);
		//CFL
		//dt=(max(dx,dy)*(max(dx,dy)*cfl))/v.max();
		dt = 0.005;
		//(Vx*dt/dx+Vy*dt/dy) <= 1
		//dt=cfl/(dx/v.max()+dy/v.max());
		//cout << "dt = " << dt << endl;
		tfinal=Nt*dt;
		//cout << "Vmax = " << v.max() << endl;

		for(int i(1); i < Ny; ++i)
			v.GetAllVX()(i,0)=0;

		Density n(Nx,Ny,L,H);
		precision rho_mer(1032.0), p_atm(1015000.0), y, g(0);

		v.GetAllVX().SetBoundaryCondition(dirichlet,0,dirichlet,0,dirichlet,1,neumann,0);
		v.GetAllVY().SetBoundaryCondition(dirichlet,0,dirichlet,0,dirichlet,0,neumann,0);
		Matrix p(Ny+1,Nx+1);
		p.SetBoundaryCondition(neumann,0,neumann,0,dirichlet,p_atm,dirichlet,p_atm);
		for(int i(0); i < Ny+1; ++i) {
			y = H-i*dy;
			for(int j(0); j < Nx+1; ++j) {
				n(i,j) = rho_mer;
				p(i,j) = p_atm + rho_mer*g*y;
			}
		}
		NavierStokes2 ns(Nx,Ny,Nt,L,H,tfinal,v,n,p);
		int nDisplay(1);
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
				ns.WriteVtk("output_test/testpoiseuille/particle" + var + var2 + to_string(i/nDisplay) + ".vtk");
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

	Matrix nule;
	Velocity v(Nx,Ny,L,H);
	v.ChampsCirculaire(L/2.0,H/2, 5.0);
	//v.ChampsUniformeVx(-1.0);
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
			n.InitialCircle(L/3.0,H/3.0,0.5,1.0);
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
				timescheme.Advance(i, tn, test1,nule);
				if((i%nDisplay)==0) {
					var=(i/nDisplay < 10 ? "0" : "");
					var2=(i/nDisplay < 100 ? "0" : "");
					/*file_out << "set terminal postscript eps enhanced color" << endl;
					file_out << "set output '" << ("animate/particle" + var + var2 + to_string(i/nDisplay) + ".eps'") << endl;
					file_out << "set pm3d map" << endl;
					file_out << ("splot 'animate/particle" + to_string(i/nDisplay) + ".dat' matrix") << endl  << endl;*/


					//n.WriteGnuPlot("animate/particle" + to_string(i/nDisplay) + ".dat");

					timescheme.GetIterate().WriteVtk("output_test/testrotateDC/testrotateDC" + var + var2 + to_string(i/nDisplay) + ".vtk", dx, dy);
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

	Matrix nule;
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
			n.InitialCircle(L/3.0,H/3.0,0.5,1.0);
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
				timescheme.Advance(i, tn, test1,nule);
				if((i%nDisplay)==0) {
					var=(i/nDisplay < 10 ? "0" : "");
					var2=(i/nDisplay < 100 ? "0" : "");
					/*file_out << "set terminal postscript eps enhanced color" << endl;
					file_out << "set output '" << ("animate/particle" + var + var2 + to_string(i/nDisplay) + ".eps'") << endl;
					file_out << "set pm3d map" << endl;
					file_out << ("splot 'animate/particle" + to_string(i/nDisplay) + ".dat' matrix") << endl  << endl;*/


					//n.WriteGnuPlot("animate/particle" + to_string(i/nDisplay) + ".dat");

					timescheme.GetIterate().WriteVtk("output_test/testboundary/dirichlet" + var + var2 + to_string(i/nDisplay) + ".vtk", dx, dy);
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



//  testrotateDC();
	//testrotatedirDC();
	testNS2();
	//testpoiseuille();
		return 0;
}
