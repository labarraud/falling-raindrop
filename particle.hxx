#ifndef FILE_PARTICLE_HXX

#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <stdint.h>
#include <string>
#include <cstring>
#include <exception>

using namespace std;


class Particle

{
private:
	Vector<Vector<double> > n;
	int Nx,Ny;
	double L,H,Delta_x,Delta_y;

public:
	Particle(int Nx,int Ny,double L,double H);
	void InitialSquare(double Xcenter,double Ycenter, double intensite);
	void WriteGnuPlot(const string& nom);
	double& Getn(int i, int j);
	Vector<Vector<double> > & Getn();
	void Setn(const Vector<Vector<double> >& n);
	void WriteVtk(const string& nom);
};

#define FILE_PARTICLE_HXX
#endif
