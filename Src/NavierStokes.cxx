/*
 * NavierStokes.cxx
 *
 *  Created on: 4 avr. 2016
 *      Author: tlariviere
 */

#include "../Include/NavierStokes.hxx"


NavierStokes::NavierStokes(const Velocity& _v, const Matrix& _rho, const Matrix& _p, int Nx, int Ny, precision _L, precision _H, precision _dt)
{
	Init(_v, _rho, _p, Nx, Ny, L, H, dt);
}

void NavierStokes::Init(const Velocity& _v, const Matrix& _rho, const Matrix& _p, int Nx, int Ny, precision _L, precision _H, precision _dt)
{
	v = _v;
	p = _p;
	rho = _rho;
	L = _L;
	H = _H;
	dx = L / Nx;
	dy = H / Ny;
	dt = _dt;
	nu = 1.007e-6; // viscosité cinématique de l'eau à 20°C
	int N(Ny+1), M(Nx+1), l(N*M);
	cond_bord_p.resize(l);
	sec_membre_p.resize(l);
	for(int i(0); i < N; ++i) {
		for(int j(1); j < M; ++j) {
			cond_bord_p[Bij(i,j,M)] = 0.0;
		}
		// Condition bords haut et bas
		cond_bord_p[Bij(i,0,M)] = 0.0;
		cond_bord_p[Bij(i,Ny,M)] = 0.0;
	}
	for(int j(0); j < M; ++j) {
		// Condition bords gauche et droite
		cond_bord_p[Bij(0,j,M)] = 0.0;
		cond_bord_p[Bij(Nx,j,M)] = 0.0;
	}
}


void NavierStokes::WriteVtk(const string& nom) const
{

}

void NavierStokes::SolveLaplacianP()
{
// construire A
	// gradient conjugué... -> Ap = b+F
	int N(Ny+1), M(Nx+1), l(N*M);
	vector<precision> div(l), x;
	Div4thOrder(v, dx, dy, N, M, div);
	for(int i(0); i < N; ++i) {
		for(int j(0); j < M; ++j) {
			sec_membre_p[Bij(i,j,M)] = rho(i,j)*div[Bij(i,j,M)]/dt + cond_bord_p[Bij(i,j,M)];
		}
	}
	GradConjLaplacian(dx, dy, N, M, 1e-6, N*M+1, x, sec_membre_p); // pression sous forme de vecteur dans x
	p.Vec2Mat(x); // on transforme x en la matrice p
}

void NavierStokes::Advance(int n, double tn)
{
	SolveLaplacianP();
}


void GradConjLaplacian(precision dx, precision dy, int N, int M, precision epsilon, int Nmax, vector<precision>& x, const vector<precision>& b)
{
	unsigned int l(N*M);
	int iter;
	vector<precision> r(l), d(l), w(l);
	precision alpha, beta, nr;

	x.resize(l);
	for(unsigned int i(0); i < l; ++i) {
		x[i] = 0.0;
	}
	Laplacian4thOrder(dx, dy, N, M, x,r);
	for(unsigned int i(0); i < l; ++i) {
		r[i] -= b[i];
	}
	d=r;
	iter = 1;
	nr = VecNorme(r);

	while(nr > epsilon && iter < Nmax) {
		Laplacian4thOrder(dx, dy, N, M, d,w);

		alpha = DotProduct(d,r)/DotProduct(d,w);

		for(unsigned int i(0); i < l; ++i) {
			x[i] -= alpha*d[i];
		}

		nr = VecNorme(r);
		beta = 1.0/(nr*nr);
		for(unsigned int i(0); i < l; ++i) {
			r[i] -= alpha*w[i];
		}
		nr = VecNorme(r);
		beta = beta*nr*nr;
		for(unsigned int i(0); i < l; ++i) {
			d[i] = r[i]+beta*d[i];
		}
	}
}

/**
 * Calcul du produit matrice vecteur : out = A.v
 * avec A la matrice du laplacien discrétisé à l'ordre 4
 */
void Laplacian4thOrder(precision dx, precision dy, int N, int M, const vector<precision>& v, vector<precision>& out)
{
	precision a(12.0*dx*dx), b(12.0*dy*dy);


	// 1er bloc
	out[Bij(0,0,M)] = (16.0*v[Bij(0,1,M)]-v[Bij(0,2,M)]-30.0*v[Bij(0,0,M)])/a
					+ (16.0*v[Bij(1,0,M)]-v[Bij(2,0,M)]-30.0*v[Bij(0,0,M)])/b;
	out[Bij(0,1,M)] = (16.0*v[Bij(0,2,M)]-v[Bij(0,3,M)]-30.0*v[Bij(0,1,M)]+16.0*v[Bij(0,0,M)])/a
					+ (16.0*v[Bij(1,1,M)]-v[Bij(2,1,M)]-30.0*v[Bij(0,1,M)])/b;
	for(int j(2); j < M-2; ++j) {
		out[Bij(0,j,M)] = (16.0*v[Bij(0,j+1,M)]-v[Bij(0,j+2,M)]-30.0*v[Bij(0,j,M)]-v[Bij(0,j-2,M)]+16.0*v[Bij(0,j-1,M)])/a
						+ (16.0*v[Bij(1,j,M)]-v[Bij(2,j,M)]-30.0*v[Bij(0,j,M)])/b;
	}
	out[Bij(0,M-2,M)] = (16.0*v[Bij(0,M-1,M)]-30.0*v[Bij(0,M-2,M)]-v[Bij(0,M-4,M)]+16.0*v[Bij(0,M-3,M)])/a
					  + (16.0*v[Bij(1,M-2,M)]-v[Bij(2,M-2,M)]-30.0*v[Bij(0,M-2,M)])/b;
	out[Bij(0,M-1,M)] = (16.0*v[Bij(0,M-2,M)]-30.0*v[Bij(0,M-1,M)]-v[Bij(0,M-3,M)])/a
					  + (16.0*v[Bij(1,M-1,M)]-v[Bij(2,M-1,M)]-30.0*v[Bij(0,M-1,M)])/b;


	// 2ème bloc
	out[Bij(1,0,M)] = (16.0*v[Bij(1,1,M)]-v[Bij(1,2,M)]-30.0*v[Bij(1,0,M)])/a
					+ (16.0*v[Bij(2,0,M)]-v[Bij(3,0,M)]-30.0*v[Bij(1,0,M)]+16.0*v[Bij(0,0,M)])/b;
	out[Bij(1,1,M)] = (16.0*v[Bij(1,2,M)]-v[Bij(1,3,M)]-30.0*v[Bij(1,1,M)]+16.0*v[Bij(1,0,M)])/a
					+ (16.0*v[Bij(2,1,M)]-v[Bij(3,1,M)]-30.0*v[Bij(1,1,M)]+16.0*v[Bij(0,1,M)])/b;
	for(int j(2); j < M-2; ++j) {
		out[Bij(1,j,M)] = (16.0*v[Bij(1,j+1,M)]-v[Bij(1,j+2,M)]-30.0*v[Bij(1,j,M)]-v[Bij(1,j-2,M)]+16.0*v[Bij(1,j-1,M)])/a
						+ (16.0*v[Bij(2,j,M)]-v[Bij(3,j,M)]-30.0*v[Bij(1,j,M)]+16.0*v[Bij(0,j,M)])/b;
	}
	out[Bij(1,M-2,M)] = (16.0*v[Bij(1,M-1,M)]-30.0*v[Bij(1,M-2,M)]-v[Bij(1,M-4,M)]+16.0*v[Bij(1,M-3,M)])/a
					  + (16.0*v[Bij(2,M-2,M)]-v[Bij(3,M-2,M)]-30.0*v[Bij(1,M-2,M)]+16.0*v[Bij(0,M-2,M)])/b;
	out[Bij(1,M-1,M)] = (16.0*v[Bij(1,M-2,M)]-30.0*v[Bij(1,M-1,M)]-v[Bij(1,M-3,M)])/a
					  + (16.0*v[Bij(2,M-1,M)]-v[Bij(3,M-1,M)]-30.0*v[Bij(1,M-1,M)]+16.0*v[Bij(0,M-1,M)])/b;


	// Milieu Matrice
	for(int i(2); i < N-2; ++i) {
		out[Bij(i,0,M)] = (16.0*v[Bij(i,1,M)]-v[Bij(i,2,M)]-30.0*v[Bij(i,0,M)])/a
						+ (16.0*v[Bij(i+1,0,M)]-v[Bij(i+2,0,M)]-30.0*v[Bij(i,0,M)]-v[Bij(i-2,0,M)]+16.0*v[Bij(i-1,0,M)])/b;
		out[Bij(i,1,M)] = (16.0*v[Bij(i,2,M)]-v[Bij(i,3,M)]-30.0*v[Bij(i,1,M)]+16.0*v[Bij(i,0,M)])/a
						+ (16.0*v[Bij(i+1,1,M)]-v[Bij(i+2,1,M)]-30.0*v[Bij(i,1,M)]-v[Bij(i-2,1,M)]+16.0*v[Bij(i-1,1,M)])/b;
		for(int j(2); j < M-2; ++j) {
			out[Bij(i,j,M)] = (16.0*v[Bij(i,j+1,M)]-v[Bij(i,j+2,M)]-30.0*v[Bij(i,j,M)]-v[Bij(i,j-2,M)]+16.0*v[Bij(i,j-1,M)])/a
							+ (16.0*v[Bij(i+1,j,M)]-v[Bij(i+2,j,M)]-30.0*v[Bij(i,j,M)]-v[Bij(i-2,j,M)]+16.0*v[Bij(i-1,j,M)])/b;
		}
		out[Bij(i,M-2,M)] = (16.0*v[Bij(i,M-1,M)]-30.0*v[Bij(i,M-2,M)]-v[Bij(i,M-4,M)]+16.0*v[Bij(i,M-3,M)])/a
						  + (16.0*v[Bij(i+1,M-2,M)]-v[Bij(i+2,M-2,M)]-30.0*v[Bij(i,M-2,M)]-v[Bij(i-2,M-2,M)]+16.0*v[Bij(i-1,M-2,M)])/b;
		out[Bij(i,M-1,M)] = (16.0*v[Bij(i,M-2,M)]-30.0*v[Bij(i,M-1,M)]-v[Bij(i,M-3,M)])/a
						  + (16.0*v[Bij(i+1,M-1,M)]-v[Bij(i+2,M-1,M)]-30.0*v[Bij(i,M-1,M)]-v[Bij(i-2,M-1,M)]+16.0*v[Bij(i-1,M-1,M)])/b;
	}


	// Avant dernier bloc
	out[Bij(N-2,0,M)] = (16.0*v[Bij(N-2,1,M)]-v[Bij(N-2,2,M)]-30.0*v[Bij(N-2,0,M)])/a
					+ (16.0*v[Bij(N-1,0,M)]-30.0*v[Bij(N-2,0,M)]-v[Bij(N-4,0,M)]+16.0*v[Bij(N-3,0,M)])/b;
	out[Bij(N-2,1,M)] = (16.0*v[Bij(N-2,2,M)]-v[Bij(N-2,3,M)]-30.0*v[Bij(N-2,1,M)]+16.0*v[Bij(N-2,0,M)])/a
					+ (16.0*v[Bij(N-1,1,M)]-30.0*v[Bij(N-2,1,M)]-v[Bij(N-4,1,M)]+16.0*v[Bij(N-3,1,M)])/b;
	for(int j(2); j < M-2; ++j) {
		out[Bij(N-2,j,M)] = (16.0*v[Bij(N-2,j+1,M)]-v[Bij(N-2,j+2,M)]-30.0*v[Bij(N-2,j,M)]-v[Bij(N-2,j-2,M)]+16.0*v[Bij(N-2,j-1,M)])/a
						+ (16.0*v[Bij(N-1,j,M)]-30.0*v[Bij(N-2,j,M)]-v[Bij(N-4,j,M)]+16.0*v[Bij(N-3,j,M)])/b;
	}
	out[Bij(N-2,M-2,M)] = (16.0*v[Bij(N-2,M-1,M)]-30.0*v[Bij(N-2,M-2,M)]-v[Bij(N-2,M-4,M)]+16.0*v[Bij(N-2,M-3,M)])/a
					  + (16.0*v[Bij(N-1,M-2,M)]-30.0*v[Bij(N-2,M-2,M)]-v[Bij(N-4,M-2,M)]+16.0*v[Bij(N-3,M-2,M)])/b;
	out[Bij(N-2,M-1,M)] = (16.0*v[Bij(N-2,M-2,M)]-30.0*v[Bij(N-2,M-1,M)]-v[Bij(N-2,M-3,M)])/a
					  + (16.0*v[Bij(N-1,M-1,M)]-30.0*v[Bij(N-2,M-1,M)]-v[Bij(N-4,M-1,M)]+16.0*v[Bij(N-3,M-1,M)])/b;


	// Dernier bloc
	out[Bij(N-1,0,M)] = (16.0*v[Bij(N-1,1,M)]-v[Bij(N-1,2,M)]-30.0*v[Bij(N-1,0,M)])/a
					+ (16.0*v[Bij(N-2,0,M)]-30.0*v[Bij(N-1,0,M)]-v[Bij(N-3,0,M)])/b;
	out[Bij(N-1,1,M)] = (16.0*v[Bij(N-1,2,M)]-v[Bij(N-1,3,M)]-30.0*v[Bij(N-1,1,M)]+16.0*v[Bij(N-1,0,M)])/a
					+ (16.0*v[Bij(N-2,1,M)]-30.0*v[Bij(N-1,1,M)]-v[Bij(N-3,1,M)])/b;
	for(int j(2); j < M-2; ++j) {
		out[Bij(N-1,j,M)] = (16.0*v[Bij(N-1,j+1,M)]-v[Bij(N-1,j+2,M)]-30.0*v[Bij(N-1,j,M)]-v[Bij(N-1,j-2,M)]+16.0*v[Bij(N-1,j-1,M)])/a
						+ (16.0*v[Bij(N-2,j,M)]-30.0*v[Bij(N-1,j,M)]-v[Bij(N-3,j,M)])/b;
	}
	out[Bij(N-1,M-2,M)] = (16.0*v[Bij(N-1,M-1,M)]-30.0*v[Bij(N-1,M-2,M)]-v[Bij(N-1,M-4,M)]+16.0*v[Bij(N-1,M-3,M)])/a
					  + (16.0*v[Bij(N-2,M-2,M)]-30.0*v[Bij(N-1,M-2,M)]-v[Bij(N-3,M-2,M)])/b;
	out[Bij(N-1,M-1,M)] = (16.0*v[Bij(N-1,M-2,M)]-30.0*v[Bij(N-1,M-1,M)]-v[Bij(N-1,M-3,M)])/a
					  + (16.0*v[Bij(N-2,M-1,M)]-30.0*v[Bij(N-1,M-1,M)]-v[Bij(N-3,M-1,M)])/b;

}

void Div4thOrder(const Velocity& v, precision dx, precision dy, int N, int M, vector<precision>& out)
{
	precision a(12.0*dx*dx), b(12.0*dy*dy);


	// 1er bloc
	out[Bij(0,0,M)] = (8.0*v.GetVX(0,1)-v.GetVX(0,2))/a
					+ (8.0*v.GetVY(1,0)-v.GetVY(2,0))/b;
	out[Bij(0,1,M)] = (8.0*v.GetVX(0,2)-v.GetVX(0,3)-8.0*v.GetVX(0,0))/a
					+ (8.0*v.GetVY(1,1)-v.GetVY(2,1))/b;
	for(int j(2); j < M-2; ++j) {
		out[Bij(0,j,M)] = (8.0*v.GetVX(0,j+1)-v.GetVX(0,j+2)+v.GetVX(0,j-2)-8.0*v.GetVX(0,j-1))/a
						+ (8.0*v.GetVY(1,j)-v.GetVY(2,j))/b;
	}
	out[Bij(0,M-2,M)] = (8.0*v.GetVX(0,M-1)+v.GetVX(0,M-4)-8.0*v.GetVX(0,M-3))/a
					+ (8.0*v.GetVY(1,M-2)-v.GetVY(2,M-2))/b;
	out[Bij(0,M-1,M)] = (v.GetVX(0,M-3)-8.0*v.GetVX(0,M-2))/a
					+ (8.0*v.GetVY(1,M-1)-v.GetVY(2,M-1))/b;


	// 2ème bloc
	out[Bij(1,0,M)] = (8.0*v.GetVX(1,1)-v.GetVX(1,2))/a
					+ (8.0*v.GetVY(2,0)-v.GetVY(3,0)-8.0*v.GetVY(0,0))/b;
	out[Bij(1,1,M)] = (8.0*v.GetVX(1,2)-v.GetVX(1,3)-8.0*v.GetVX(1,0))/a
					+ (8.0*v.GetVY(2,1)-v.GetVY(3,1)-8.0*v.GetVY(0,1))/b;
	for(int j(2); j < M-2; ++j) {
		out[Bij(1,j,M)] = (8.0*v.GetVX(1,j+1)-v.GetVX(1,j+2)+v.GetVX(1,j-2)-8.0*v.GetVX(1,j-1))/a
						+ (8.0*v.GetVY(2,j)-v.GetVY(3,j)-8.0*v.GetVY(0,j))/b;
	}
	out[Bij(1,M-2,M)] = (8.0*v.GetVX(1,M-1)+v.GetVX(1,M-4)-8.0*v.GetVX(1,M-3))/a
					+ (8.0*v.GetVY(2,M-2)-v.GetVY(3,M-2)-8.0*v.GetVY(0,M-2))/b;
	out[Bij(1,M-1,M)] = (v.GetVX(1,M-3)-8.0*v.GetVX(1,M-2))/a
					+ (8.0*v.GetVY(2,M-1)-v.GetVY(3,M-1)-8.0*v.GetVY(0,M-1))/b;


	// Milieu Matrice
	for(int i(2); i < N-2; ++i) {
		out[Bij(i,0,M)] = (8.0*v.GetVX(i,1)-v.GetVX(i,2))/a
						+ (8.0*v.GetVY(i+1,0)-v.GetVY(i+2,0)+v.GetVY(i-2,0)-8.0*v.GetVY(i-1,0))/b;
		out[Bij(i,1,M)] = (8.0*v.GetVX(i,2)-v.GetVX(i,3)-8.0*v.GetVX(i,0))/a
						+ (8.0*v.GetVY(i+1,1)-v.GetVY(i+2,1)+v.GetVY(i-2,1)-8.0*v.GetVY(i-1,1))/b;
		for(int j(2); j < M-2; ++j) {
			out[Bij(i,j,M)] = (8.0*v.GetVX(i,j+1)-v.GetVX(i,j+2)+v.GetVX(i,j-2)-8.0*v.GetVX(i,j-1))/a
							+ (8.0*v.GetVY(i+1,j)-v.GetVY(i+2,j)+v.GetVY(i-2,j)-8.0*v.GetVY(i-1,j))/b;
		}
		out[Bij(i,M-2,M)] = (8.0*v.GetVX(i,M-1)+v.GetVX(i,M-4)-8.0*v.GetVX(i,M-3))/a
						+ (8.0*v.GetVY(i+1,M-2)-v.GetVY(i+2,M-2)+v.GetVY(i-2,M-2)-8.0*v.GetVY(i-1,M-2))/b;
		out[Bij(i,M-1,M)] = (v.GetVX(i,M-3)-8.0*v.GetVX(i,M-2))/a
						+ (8.0*v.GetVY(i+1,M-1)-v.GetVY(i+2,M-1)+v.GetVY(i-2,M-1)-8.0*v.GetVY(i-1,M-1))/b;
	}


	// Avant dernier bloc
	out[Bij(N-2,0,M)] = (8.0*v.GetVX(N-2,1)-v.GetVX(N-2,2))/a
					+ (8.0*v.GetVY(N-1,0)+v.GetVY(N-4,0)-8.0*v.GetVY(N-3,0))/b;
	out[Bij(N-2,1,M)] = (8.0*v.GetVX(N-2,2)-v.GetVX(N-2,3)-8.0*v.GetVX(N-2,0))/a
					+ (8.0*v.GetVY(N-1,1)+v.GetVY(N-4,1)-8.0*v.GetVY(N-3,1))/b;
	for(int j(2); j < M-2; ++j) {
		out[Bij(N-2,j,M)] = (8.0*v.GetVX(N-2,j+1)-v.GetVX(N-2,j+2)+v.GetVX(N-2,j-2)-8.0*v.GetVX(N-2,j-1))/a
						+ (8.0*v.GetVY(N-1,j)+v.GetVY(N-4,j)-8.0*v.GetVY(N-3,j))/b;
	}
	out[Bij(N-2,M-2,M)] = (8.0*v.GetVX(N-2,M-1)+v.GetVX(N-2,M-4)-8.0*v.GetVX(N-2,M-3))/a
					+ (8.0*v.GetVY(N-1,M-2)+v.GetVY(N-4,M-2)-8.0*v.GetVY(N-3,M-2))/b;
	out[Bij(N-2,M-1,M)] = (v.GetVX(N-2,M-3)-8.0*v.GetVX(N-2,M-2))/a
					+ (8.0*v.GetVY(N-1,M-1)+v.GetVY(N-4,M-1)-8.0*v.GetVY(N-3,M-1))/b;


	// Dernier bloc
	out[Bij(N-1,0,M)] = (8.0*v.GetVX(N-1,1)-v.GetVX(N-1,2))/a
					+ (+v.GetVY(N-3,0)-8.0*v.GetVY(N-2,0))/b;
	out[Bij(N-1,1,M)] = (8.0*v.GetVX(N-1,2)-v.GetVX(N-1,3)-8.0*v.GetVX(N-1,0))/a
					+ (v.GetVY(N-3,1)-8.0*v.GetVY(N-2,1))/b;
	for(int j(2); j < M-2; ++j) {
		out[Bij(N-1,j,M)] = (8.0*v.GetVX(N-1,j+1)-v.GetVX(N-1,j+2)+v.GetVX(N-1,j-2)-8.0*v.GetVX(N-1,j-1))/a
						+ (v.GetVY(N-3,j)-8.0*v.GetVY(N-2,j))/b;
	}
	out[Bij(N-1,M-2,M)] = (8.0*v.GetVX(N-1,M-1)+v.GetVX(N-1,M-4)-8.0*v.GetVX(N-1,M-3))/a
					+ (v.GetVY(N-3,M-2)-8.0*v.GetVY(N-2,M-2))/b;
	out[Bij(N-1,M-1,M)] = (v.GetVX(N-1,M-3)-8.0*v.GetVX(N-1,M-2))/a
					+ (v.GetVY(N-3,M-1)-8.0*v.GetVY(N-2,M-1))/b;

}


