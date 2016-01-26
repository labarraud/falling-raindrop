#ifndef FILE_VELOCITY_CXX

#include "particle.hxx"

//	Vector<Vector<double> > n;
//	int Nx,Ny;
//	double L,H,Delta_x,Delta_y;

namespace linalg
{

inline Velocity::Velocity(int Nx,int Ny,double L,double H)
	{

		this->Nx=Nx;
		this->Ny=Ny;

		this->n.Reallocate(Nx+1);

		for (int i=0; i<Nx+1;i++)
			this->n(i).Reallocate(Ny+1);

		for (int i=0; i<Nx+1;i++)
		{
			for (int i=0; i<Nx+1;i++)
				(this->n)(i)(j);
		}

		this->L=L;
		this->H=H;
		this->Delta_x=L/Nx;
		this->Delta_y=H/Ny;


	}

inline void Velocity::InitialSquare(double Xcenter,double Ycenter, double intensite)
	{
		Vector<double> v(2);
		for (int i=0; i<Nx+1;i++)
		{
			for (int j=0; j<Ny+1;j++)
			{
				v(1) = (j * Delta_y)- Ycenter;
				v(0) = (i * Delta_x)- Xcenter;

				if()

			}
		}
	}




inline void Velocity::WriteGnuPlot(const string& nom)
	{
			ofstream file_out(nom.data());
		  file_out.precision(15);
		  double x,y;
		  for (int i = 0; i < Nx+1; i++)
				for (int j=0; j<Ny+1;j++)
				{
					x=i*this->Delta_x;
					y=j*this->Delta_y;
					file_out << x << " "<< y << " " << VX(i)(j) << " " << VY(i)(j) << '\n';
				}
		  file_out.close();
	}


inline	double& Velocity::GetVX(int i, int j)
	{
		return VX(i)(j);
	}

inline	double& Velocity::GetVY(int i,int j)
	{
		return VY(i)(j);
	}
}

#define FILE_VELOCITY_CXX
#endif
