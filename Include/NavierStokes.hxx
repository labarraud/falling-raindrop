/*
 * NavierStokes.hxx
 *
 *  Created on: 4 avr. 2016
 *      Author: tlariviere
 */

#ifndef INCLUDE_NAVIERSTOKES_HXX_
#define INCLUDE_NAVIERSTOKES_HXX_

#include "Matrix.hxx"
#include "velocity.hxx"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class NavierStokes
{
public:
	NavierStokes(const Velocity& _v, const Matrix& _rho, const Matrix& _p, int Nx, int Ny, precision _L, precision _H, precision _dt);
	void Init(const Velocity& _v, const Matrix& _rho, const Matrix& _p, int Nx, int Ny, precision _L, precision _H, precision _dt);
	void WriteVtk(const string& nom) const;
	void SolveLaplacianP();
	void Advance(int n, double tn);
private:
	Velocity v;
	Matrix rho; // masse volumique
	Matrix p; // pression
	vector<precision> cond_bord_p;
	vector<precision> sec_membre_p;
	int Nx, Ny;
	precision L,H,dx,dy,dt,nu;
};

void GradConjLaplacian(precision dx, precision dy, int N, int M, precision epsilon, int Nmax, vector<precision>& x, const vector<precision>& b);
void Laplacian4thOrder(precision dx, precision dy, int N, int M, const vector<precision>& v, vector<precision>& out);
void Div4thOrder(const Velocity& v, precision dx, precision dy, int N, int M, vector<precision>& out);


#endif /* INCLUDE_NAVIERSTOKES_HXX_ */
