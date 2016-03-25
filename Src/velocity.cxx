#ifndef FILE_VELOCITY_CXX

#include "../Include/velocity.hxx"

//	Vector<Vector<double>> VX,VY;
	//int Nx,Ny;
	//double L,H,Delta_x,Delta_y;


Velocity::Velocity(){}
Velocity::Velocity(int Nx,int Ny,double L,double H)
	{

		this->Nx=Nx;
		this->Ny=Ny;

		this->VX.Reallocate(Nx+1,Ny+1);
		this->VY.Reallocate(Nx+1,Ny+1);

		this->L=L;
		this->H=H;
		this->Delta_x=L/Nx;
		this->Delta_y=H/Ny;


	}

void Velocity::ChampsCirculaire(double Xcenter,double Ycenter, double intensite)
	{
		vector<precision> v(2);
		for (int i=0; i<Ny+1;i++)
		{
			for (int j=0; j<Nx+1;j++)
			{
				v[1] = (i * Delta_y)- Ycenter;
				v[0] = (j * Delta_x)- Xcenter;
				VX(i,j) = intensite*v[1];
				VY(i,j) = -intensite*v[0];
			}
		}
	}




void Velocity::WriteGnuPlot(const string& nom)
	{
			ofstream file_out(nom.data());
		  file_out.precision(15);
		  double x,y;
		  for (int i = 0; i < Ny+1; i++)
				for (int j=0; j<Nx+1;j++)
				{
					x=j*this->Delta_x;
					y=i*this->Delta_y;
					file_out << x << " "<< y << " " << VX(i,j) << " " << VY(i,j) << '\n';
				}
		  file_out.close();
	}


precision& Velocity::GetVX(int i, int j)
	{
		return VX(i,j);
	}

precision& Velocity::GetVY(int i,int j)
	{
		return VY(i,j);
	}


precision Velocity::max()
{
	double max= abs(VX(0,0));
	for(int i =0; i<Nx ; i++)
	{	
		for(int j =0; j<Ny ; j++)
		{
			if(abs(VX(i,j)) > max) {
				max = abs(VX(i,j));
			}
			if(abs(VY(i,j)) > max) {
				max = abs(VY(i,j));
			}
		}	
	}
	return max;
}

void Velocity::ChampsUniformeVx(double intensite)
{
			for (int i=0; i<Ny+1;i++)
			{
				for (int j=0; j<Nx+1;j++)
				{
					VX(i,j) = intensite;
					VY(i,j) = 0;
				}
			}

}

void Velocity::ChampsUniforme(double intensite)
{
			for (int i=0; i<Ny+1;i++)
			{
				for (int j=0; j<Nx+1;j++)
				{
					VX(i,j) = intensite;
					VY(i,j) = intensite;
				}
			}

}

#define FILE_VELOCITY_CXX
#endif
