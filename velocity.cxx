#ifndef FILE_VELOCITY_CXX

#include "velocity.hxx"


//	Vector<Vector<double>> VX,VY;
	//int Nx,Ny;
	//double L,H,Delta_x,Delta_y;

	Velocity::Velocity(int Nx,int Ny,double L,double H)
	{
		this.VX.realocate(Nx+1);
		this.VX.realocate(Nx+1);

		for (int i=0; i<Nx+1;i++)
		{
			this.VX(i).realocate(Ny+1)
			this.VY(i).realocate(Ny+1)
		}

		this.Delta_x=L/Nx;
		this.Delta_y=H/Ny;


	}

	void Velocity::ChampsCirculaire(double Xcenter,double Ycenter, double intensite)
	{
		Vector(double) v(2);
		for (int i=0; i<Nx+1;i++)
		{
			for (int j=0; j<Ny+1;j++)
			{
				v(1) = Xcenter - (i * Delta_y);
				v(0) = Ycenter - (j * Delta_x);
				VX(i)(j) = intensite*v(1);
				VY(i)(j) = -intensite*v(0);
			}
		}
	}




	void Velocity::WriteGnuPlot(const string& nom)
	{
		  ofstream file_out(nom.data());
		  file_out.precision(15);
		  double x,y;
		  for (int i = 0; i < Nx; i++)
				for (int j=0; j<Ny+1;j++)
				{
					x=i*Delta_x;
					y=j*Delta_y;
					file_out << x << " "<< y << " " << VX(i)(j) << " " << VY(i)(j) << '\n';
				}
		  file_out.close();
	}


	double& Velocity::GetVX(int i, j)
	{
		return VX(i,j);
	}

	double& Velocity::GetVY(int i, j)
	{
		return VY(i,j);
	}


#define FILE_VELOCITY_CXX
#endif
