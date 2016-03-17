#ifndef FILE_PARTICLE_CXX

#include "particle.hxx"

//	Vector<Vector<double> > n;
//	int Nx,Ny;
//	double L,H,Delta_x,Delta_y;


Particle::Particle(int Nx,int Ny,double L,double H)
	{

		this->Nx=Nx;
		this->Ny=Ny;

		this->n.Reallocate(Nx+1,Ny+1);



		for (int i=0; i<Nx+1;i++)
		{
			for (int j=0; j<Ny+1;j++)
				n(i,j)=0;
		}

		this->L=L;
		this->H=H;
		this->Delta_x=L/Nx;
		this->Delta_y=H/Ny;
	}

void Particle::InitialSquare(precision Xcenter,precision Ycenter, precision intensite)
	{
		vector<precision> v(2);
		for (int i=0; i<Ny+1;i++)
		{
			for (int j=0; j<Nx+1;j++)
			{
				v[1] = std::abs(Ycenter - (i * Delta_y));
				v[0] = std::abs(Xcenter-(j * Delta_x));

				if(std::max(v[0],v[1])<(intensite*intensite))
						n(i,j)=1;

			}
		}
	}




void Particle::WriteGnuPlot(const string& nom)
	{
		ofstream file_out(nom.data());
		file_out.precision(15);
	//	double x,y;
		for (int i = 0; i < Ny+1; i++)
			{
			for (int j=0; j<Nx+1;j++)
			{
				file_out << n(i,j)<< " ";
			}
				file_out << "\n";
		  }
		  file_out.close();
	}



Matrix & Particle::Getn()
{
	return n;
}


void Particle::Setn(const Matrix& n)
{
	this->n=n;
	
}

void Particle::WriteVtk(const string& nom)
	{
		  ofstream file_out(nom.data());
		  file_out << "# vtk DataFile Version 3.1\n";
		  file_out << "2-D mesh\n";
		  file_out << "ASCII\n\n";
		  file_out << "DATASET STRUCTURED_GRID\n";
		  file_out << "DIMENSIONS " << (Nx+1) << " " << (Ny+1) << " 1\n";
		  file_out << "SPACING " << Delta_x << " " << Delta_y << " 1.0\n\n";
		  file_out << "ORIGIN 0.0 0.0 0.0\n";
		  file_out << "POINT_DATA " << ((Nx+1)*(Ny+1)) << " \n";
		  file_out << "SCALARS Concentration FLOAT 1\n";
		  file_out << "LOOKUP_TABLE default\n";
		  for (int i = 0; i < Ny+1; i++) {
				for (int j=0; j<Nx+1;j++) {
//					file_out << n(i)(j) << '\n';
				}
		  }
		  file_out.close();
	}

#define FILE_PARTICLE_CXX
#endif
