#ifndef FILE_PARTICLE_HXX

#include"Matrix.hxx"
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
	Matrix n;
	int Nx,Ny;
	precision L,H,Delta_x,Delta_y;

public:
	Particle(int Nx,int Ny,precision L,precision H);
	void InitialSquare(precision Xcenter,precision Ycenter, precision intensite);
	void WriteGnuPlot(const string& nom);
	double& Getn(int i, int j);
	Matrix& Getn();
	void Setn(const Matrix& n);
	void WriteVtk(const string& nom);
};

#define FILE_PARTICLE_HXX
#endif
