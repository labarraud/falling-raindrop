#ifndef FILE_VELOCITY_CXX

#include "velocity.hxx"

#include <fstream>
//	Vector<Vector<double>> VX,VY;
	//int Nx,Ny;
	//double L,H,Delta_x,Delta_y;

namespace linalg
{

inline Velocity::Velocity(int Nx,int Ny,double L,double H)
	{
		this->VX.Reallocate(Nx+1);
		this->VY.Reallocate(Ny+1);

		for (int i=0; i<Nx+1;i++)
		{
			this->VX(i).Reallocate(Nx+1);
			this->VY(i).Reallocate(Ny+1);
		}

		this->Delta_x=L/Nx;
		this->Delta_y=H/Ny;


	}

inline void Velocity::ChampsCirculaire(double Xcenter,double Ycenter, double intensite)
	{
		Vector<double> v(2);
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




inline void Velocity::WriteGnuPlot(const string& nom)
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
