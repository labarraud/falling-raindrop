#ifndef FILE_VELOCITY_HXX

#include"define.hxx"
#include"Matrix.hxx"

#include <vector>
#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <stdint.h>
#include <string>
#include <cstring>
#include <exception>

using namespace std;


class Velocity

{
private:
	Matrix VX,VY;
	int Nx,Ny;
	precision L,H,Delta_x,Delta_y;

public:
	Velocity();
	Velocity(int Nx,int Ny,precision L,precision H);
	void ChampsCirculaire(precision Xcenter,precision Ycenter, precision intensite);
	void ChampsCircle(precision Xcenter,precision Ycenter, precision radius, precision vx, precision vy);
	void ChampsUniformeVx(precision intensite);
	void ChampsUniformeVy(precision intensite);
	void ChampsUniforme(precision intensite);
	void WriteGnuPlot(const string& nom);
	Matrix& GetAllVX();
	Matrix& GetAllVY();
	void SetAllVX(const Matrix& _VX);
	void SetAllVY(const Matrix& _VY);
	precision& GetVX(int i, int j);
	precision& GetVY(int i, int j);
	precision GetVX(int i, int j) const;
	precision GetVY(int i, int j) const;
	precision max();

};

#define FILE_VELOCITY_HXX
#endif
