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
#include "DiffusionConvectionProblem.hxx"
#include "TimeScheme.hxx"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class NavierStokes : public VirtualOdeSystem
{
public:
	NavierStokes();
	NavierStokes(int _Nx,int _Ny,int _Nt,precision _L,precision _H,precision tfinal,Velocity& _v,Density& _rho,const Matrix& _p);
	void SetPressure(const Matrix& _p);
	virtual void SetInitialCondition(int _Nx,int _Ny,int _Nt,precision _L,precision _H,precision tfinal,Velocity& _v,Density& _rho);
	void WriteVtk(const string& nom) const;
	void SolveLaplacianP();
	void Advance(int n, double tn);
	virtual void AddFunction(precision alpha, const Matrix& rho, precision t, Matrix& y);
	precision UpwindY(precision dt, precision a, int i, int j, const Matrix& u);
	precision SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u);
	precision uij(int i, int j, const Matrix& u);
	precision p_bord_droit(int i,int j) const;
	precision p_bord_gauche(int i,int j) const;
	precision p_bord_haut(int i,int j) const;
	precision p_bord_bas(int i,int j) const;
	precision vx_bord_droit(int i,int j, const Matrix& u) const;
	precision vx_bord_gauche(int i,int j, const Matrix& u) const;
	precision vx_bord_haut(int i,int j, const Matrix& u) const;
	precision vx_bord_bas(int i,int j, const Matrix& u) const;
	precision vy_bord_droit(int i,int j, const Matrix& u) const;
	precision vy_bord_gauche(int i,int j, const Matrix& u) const;
	precision vy_bord_haut(int i,int j, const Matrix& u) const;
	precision vy_bord_bas(int i,int j, const Matrix& u) const;
private:
	Velocity v;
	Density rho; // masse volumique
	vector<precision> p; // pression
	vector<precision> cond_bord_p;
	vector<precision> cond_bord_gradx_p;
	vector<precision> cond_bord_grady_p;
	vector<precision> cond_bord_div_v;
	vector<precision> cond_bord_v;
	vector<precision> sec_membre_p;
	vector<precision> sec_membre_vx, sec_membre_vy;
	int Nx, Ny, Nt;
	precision L,H,dx,dy,dt,nu,g;
	LowStorageRungeKuttaIterator timescheme_x, timescheme_y;
	bool vx;
};

void GradConjLaplacian(precision dx, precision dy, int N, int M, precision epsilon, int Nmax, vector<precision>& x, const vector<precision>& b);
void Laplacian4thOrder(precision dx, precision dy, int N, int M, const vector<precision>& v, vector<precision>& out);
void Div4thOrder(const Velocity& v, precision dx, precision dy, int N, int M, vector<precision>& out);
void Gradx4thOrder(precision dx, int N, int M, const vector<precision>& v, vector<precision>& out);
void Grady4thOrder(precision dy, int N, int M, const vector<precision>& v, vector<precision>& out);


#endif /* INCLUDE_NAVIERSTOKES_HXX_ */
