#ifndef FILE_VELOCITY_CXX

#include "velocity.hxx"

//	Vector<Vector<double>> VX,VY;
	//int Nx,Ny;
	//double L,H,Delta_x,Delta_y;



inline Velocity::Velocity(int Nx,int Ny,double L,double H)
	{

		this->Nx=Nx;
		this->Ny=Ny;

		this->VX.Reallocate(Nx+1);
		this->VY.Reallocate(Nx+1);

		for (int i=0; i<Nx+1;i++)
		{
			this->VX(i).Reallocate(Ny+1);
			this->VY(i).Reallocate(Ny+1);
		}

		this->L=L;
		this->H=H;
		this->Delta_x=L/Nx;
		this->Delta_y=H/Ny;


	}

inline void Velocity::ChampsCirculaire(double Xcenter,double Ycenter, double intensite)
	{
		Vector<double> v(2);
		for (int i=0; i<Ny+1;i++)
		{
			for (int j=0; j<Nx+1;j++)
			{
				v(1) = (i * Delta_y)- Ycenter;
				v(0) = (j * Delta_x)- Xcenter;
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
		  for (int i = 0; i < Ny+1; i++)
				for (int j=0; j<Nx+1;j++)
				{
					x=j*this->Delta_x;
					y=i*this->Delta_y;
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


inline double Velocity::max()
{
	double max= VX(0)(0);
	for(int i =0; i<Nx ; i++)
	{	
		for(int j =0; j<Ny ; j++)
		{
			if(VX(i)(j) > max) {
				max = VX(i)(j);
			}
			if(VY(i)(j) > max) {
				max = VY(i)(j);
			}
		}	
	}
	return max;
}

inline void Velocity::ChampsUniformeVx(double intensite)
{
			for (int i=0; i<Ny+1;i++)
			{
				for (int j=0; j<Nx+1;j++)
				{
					VX(i)(j) = intensite;
					VY(i)(j) = 0;
				}
			}

}

inline void Velocity::ChampsUniforme(double intensite)
{
			for (int i=0; i<Ny+1;i++)
			{
				for (int j=0; j<Nx+1;j++)
				{
					VX(i)(j) = intensite;
					VY(i)(j) = intensite;
				}
			}

}

#define FILE_VELOCITY_CXX
#endif
