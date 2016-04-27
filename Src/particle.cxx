#ifndef FILE_PARTICLE_CXX

#include "../Include/particle.hxx"

//	Vector<Vector<double> > n;
//	int M,N;
//	double L,H,Delta_x,Delta_y;

Density::Density()
	:	Matrix()
{}

Density::Density(int Nx,int Ny,double L,double H)
	:	Matrix(Ny, Nx)
	{
		this->L=L;
		this->H=H;
		this->Delta_x=L/Nx;
		this->Delta_y=H/Ny;
	}

void Density::InitialSquare(precision Xcenter,precision Ycenter, precision intensite)
	{
		vector<precision> v(2);
		for (int i=0; i<N;i++)
		{
			for (int j=0; j<M;j++)
			{
				v[1] = std::abs(Ycenter - (i * Delta_y));
				v[0] = std::abs(Xcenter-(j * Delta_x));

				if(std::max(v[0],v[1])<(intensite*intensite))
					(*this)(i,j)=1;

			}
		}
	}

void Density::InitialGauss(precision Xcenter,precision Ycenter, precision intensite)
	{
		vector<precision> v(2);
		for (int i=0; i<N;i++)
		{
			for (int j=0; j<M;j++)
			{

				if(std::max(v[0],v[1])<(intensite*intensite))
					(*this)(i,j)=std::exp(-(Ycenter - (i * Delta_y))*(Ycenter - (i * Delta_y))-(Xcenter-(j * Delta_x))*(Xcenter-(j * Delta_x)));

			}
		}
	}

void Density::InitialCircle(precision Xcenter,precision Ycenter, precision intensite)
	{
		precision var;
		for (int i=0; i<N;i++)
		{
			for (int j=0; j<M;j++)
			{
				var=(Xcenter-(j * Delta_x))*(Xcenter-(j * Delta_x))+(Ycenter - (i * Delta_y))*(Ycenter - (i * Delta_y));

				if(var<(intensite*intensite))
					(*this)(i,j)=1;

			}
		}
	}

#define FILE_PARTICLE_CXX
#endif
