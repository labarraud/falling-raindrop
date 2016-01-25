#ifndef FILE_VELOCITY_CXX




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

		this.Delta_x=L/(Nx+1);
		this.Delta_y=H/(Ny+1);


	}

	void Velocity::ChampsCirculaire(double Xcenter,double Ycenter, double intensite)
	{
		for (int i=0; i<Nx+1;i++)
		{
			for (int j=0; j<Ny+1;j++)
			{
				this.Vx(i)(j)=;
			}
		}
	}




	void WriteGnuPlot(const string& nom);

};

#define FILE_VELOCITY_CXX
#endif
