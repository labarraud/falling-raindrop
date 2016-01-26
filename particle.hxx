#ifndef FILE_VELOCITY_HXX

#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <stdint.h>
#include <string>
#include <cstring>
#include <exception>

using namespace std;

namespace linalg
{
class Particule

{
private:
	Vector<Vector<double> > n;
	int Nx,Ny;
	double L,H,Delta_x,Delta_y;

public:
	Particule(int Nx,int Ny,double L,double H);
	void InitialSquare(double Xcenter,double Ycenter, double intensite);
	void WriteGnuPlot(const string& nom);
	double& Getn(int i, int j);

};
}
#define FILE_VELOCITY_HXX
#endif
