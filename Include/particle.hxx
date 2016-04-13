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


class Density : public Matrix
{
private:
	precision L,H,Delta_x,Delta_y;

public:
	Density();
	Density(int Nx,int Ny,precision L,precision H);
	void InitialSquare(precision Xcenter,precision Ycenter, precision intensite);
	void InitialCircle(precision Xcenter,precision Ycenter, precision intensite);
	void InitialGauss(precision Xcenter,precision Ycenter, precision intensite);
};

#define FILE_PARTICLE_HXX
#endif
