/*
 * NavierStokes2.hxx
 *
 *  Created on: 26 avr. 2016
 *      Author: tlariviere
 */

#ifndef INCLUDE_NAVIERSTOKES2_HXX_
#define INCLUDE_NAVIERSTOKES2_HXX_

#include "Matrix.hxx"
#include "velocity.hxx"
#include "DiffusionConvectionProblem.hxx"
#include "spacescheme.hxx"
#include "TimeScheme.hxx"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

class NavierStokes2
{
public:
	NavierStokes2();
	NavierStokes2(int _Nx,int _Ny,int _Nt,precision _L,precision _H,precision tfinal,Velocity& _v,Density& _rho,const Matrix& _p);
	void SetPressure(const Matrix& _p);
	virtual void SetInitialCondition(int _Nx,int _Ny,int _Nt,precision _L,precision _H,precision tfinal,Velocity& _v,Density& _rho);
	void WriteVtk(const string& nom) const;
	void SolveLaplacianP();
	void Advance(int n, double tn);

private:
	bool firstIter;
	Velocity v;
	Density rho; // masse volumique
	Matrix p; // pression
	Matrix cond_bord_p;
	Matrix sec_membre_p;
	Matrix sec_membre_vx, sec_membre_vy;
	int Nx, Ny, Nt;
	precision L,H,dx,dy,dt,nu,g,rho_mer,p_atm;
	LowStorageRungeKuttaIterator timescheme_x, timescheme_y;
	UpwindDCOrder2 spacescheme;
};

void GradConjLaplacian2(precision dx, precision dy, int N, int M, precision epsilon, int Nmax, Matrix& x, Matrix& b);
void Laplacian2thOrder(precision dx, precision dy, int N, int M, const Matrix& v, Matrix& out);
void Div2thOrder(const Velocity& v, precision dx, precision dy, int N, int M, Matrix& out);
void Gradx2thOrder(precision dx, int N, int M, const Matrix& v, Matrix& out);
void Grady2thOrder(precision dy, int N, int M, const Matrix& v, Matrix& out);


#endif /* INCLUDE_NAVIERSTOKES2_HXX_ */
