#ifndef FILE_PARTICLE_CXX

#include "particle.hxx"

//	Vector<Vector<double> > n;
//	int Nx,Ny;
//	double L,H,Delta_x,Delta_y;


inline Particle::Particle(int Nx,int Ny,double L,double H)
	{

		this->Nx=Nx;
		this->Ny=Ny;

		this->n.Reallocate(Nx+1);

		for (int i=0; i<Nx+1;i++)
			this->n(i).Reallocate(Ny+1);


		for (int i=0; i<Nx+1;i++)
		{
			for (int j=0; j<Ny+1;j++)
				n(i)(j)=0;
		}

		this->L=L;
		this->H=H;
		this->Delta_x=L/Nx;
		this->Delta_y=H/Ny;


	}

inline void Particle::InitialSquare(double Xcenter,double Ycenter, double intensite)
	{
		Vector<double> v(2);
		for (int i=0; i<Ny+1;i++)
		{
			for (int j=0; j<Nx+1;j++)
			{
				v(1) = Ycenter - (i * Delta_y);
				v(0) = Xcenter-(j * Delta_x);

				if((v(0)*v(0)+v(1)*v(1))<(intensite*intensite))
						n(i)(j)=1;

			}
		}
	}




inline void Particle::WriteGnuPlot(const string& nom)
	{
			ofstream file_out(nom.data());
		  file_out.precision(15);
	//	  double x,y;
		  for (int i = 0; i < Ny+1; i++)
		  {
				for (int j=0; j<Nx+1;j++)
				{
					file_out << n(i)(j)<< " ";
				}
				file_out << "\n";
		  }
		  file_out.close();
	}



inline 	Vector<Vector<double> > & Particle::Getn()
{
	return n;
}


inline 	void Particle::Setn(const Vector<Vector<double> >& n)
{
	this->n=n;
	
}

#define FILE_PARTICLE_CXX
#endif
