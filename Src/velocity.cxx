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

		this->VX.Reallocate(Ny+1,Nx+1);
		this->VY.Reallocate(Ny+1,Nx+1);

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


void Velocity::ChampsCircle(precision Xcenter,precision Ycenter, precision radius, precision vx, precision vy)
	{
		precision var;
		for (int i=0; i<Ny+1;i++)
		{
			for (int j=0; j<Nx+1;j++)
			{
				var=(Xcenter-(j * Delta_x))*(Xcenter-(j * Delta_x))+(Ycenter - (i * Delta_y))*(Ycenter - (i * Delta_y));

				if(var<(radius*radius)) {
					VX(i,j) = vx;
					VY(i,j) = vy;
				} else {
					VX(i,j) = 0.0;
					VY(i,j) = 0.0;
				}

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


Matrix& Velocity::GetAllVX()
{
	return VX;
}

Matrix& Velocity::GetAllVY()
{
	return VY;
}

void Velocity::SetAllVX(const Matrix& _VX)
{
	VX = _VX;
}

void Velocity::SetAllVY(const Matrix& _VY)
{
	VY = _VY;
}

precision& Velocity::GetVX(int i, int j)
	{
		return VX(i,j);
	}

precision& Velocity::GetVY(int i,int j)
	{
		return VY(i,j);
	}


precision Velocity::GetVX(int i, int j) const
	{
		return VX(i,j);
	}

precision Velocity::GetVY(int i,int j) const
	{
		return VY(i,j);
	}


precision Velocity::max()
{
	double max= abs(VX(0,0));
	for(int i =0; i<Ny+1 ; i++)
	{	
		for(int j =0; j<Nx+1 ; j++)
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

void Velocity::ChampsUniformeVy(double intensite)
{
			for (int i=0; i<Ny+1;i++)
			{
				for (int j=0; j<Nx+1;j++)
				{
					VX(i,j) = 0;
					VY(i,j) = intensite;
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
