/*
 * NavierStokes.cxx
 *
 *  Created on: 4 avr. 2016
 *      Author: tlariviere
 */

#include "../Include/NavierStokes2.hxx"

NavierStokes2::NavierStokes2()
{ }

NavierStokes2::NavierStokes2(int _Nx,int _Ny,int _Nt,precision _L,precision _H,precision tfinal,Velocity& _v,Density& _rho,const Matrix& _p)
{
	SetInitialCondition(_Nx,_Ny,_Nt,_L,_H,tfinal,_v,_rho);
	p=_p;
}

void NavierStokes2::SetInitialCondition(int _Nx,int _Ny,int _Nt,precision _L,precision _H,precision tfinal,Velocity& _v,Density& _rho)
{
	firstIter=true;
	v = _v;
	rho = _rho;
	L = _L;
	H = _H;
	Nx = _Nx;
	Ny = _Ny;
	Nt = _Nt;
	dx = L / _Nx;
	dy = H / _Ny;
	dt = tfinal/Nt;
	//nu = 1.007e-6; // viscosité cinématique de l'eau douce à 20°C
	nu = 2.8e-6; // viscosité cinématique de l'eau salée à 20°C, 50 kg/m^3 de NaCl
	rho_douce = 1000.0;
	D_nacl = 2.0e-7;
	rho_mer = 1032.0;
	p_atm = 1015000.0;
	g = 9.81;
	int N(Ny+1), M(Nx+1);
	cond_bord_p.Reallocate(N,M);
	sec_membre_p.Reallocate(N,M);
	sec_membre_vx.Reallocate(N,M);
	sec_membre_vy.Reallocate(N,M);
	zero.Reallocate(N,M);
	zero.Zero();
	spacescheme.SetInitialCondition(Nx,Ny,Nt,L,H,tfinal,v,rho);
	timescheme_x.SetInitialCondition(0,dt,v.GetAllVX(),spacescheme);
	timescheme_y.SetInitialCondition(0,dt,v.GetAllVY(),spacescheme);
	timescheme_n.SetInitialCondition(0,dt,rho,spacescheme);
}

void NavierStokes2::SetPressure(const Matrix& _p)
{
	p=_p;
}


void NavierStokes2::WriteVtk(const string& nom) const
{
	int N(Ny+1), M(Nx+1);
	ofstream file_out(nom.data());
	file_out << "# vtk DataFile Version 2.0\n";
	file_out << "Titre\n";
	file_out << "ASCII\n";
	file_out << "DATASET STRUCTURED_POINTS\n";
	file_out << "DIMENSIONS " << M << " " << N << " 1\n";
	file_out << "ORIGIN 0.0 0.0 0.0\n";
	file_out << "SPACING " << dx << " " << dy << " 0.0\n";
	file_out << "POINT_DATA " << (M*N) << " \n";
	file_out << "VECTORS velocity float\n";
	for (int i = 0; i < N; i++) {
		for (int j=0; j<M;j++) {
			if((j%5)==0) {
				file_out << '\n';
			}
			file_out << v.GetVX(i,j) << ' ' << v.GetVY(i,j) << " 0\t";
		}
	}
	file_out << "\nSCALARS pressure float\n";
	file_out << "LOOKUP_TABLE default";
	for (int i = 0; i < N; i++) {
		for (int j=0; j<M;j++) {
			if((j%10)==0) {
				file_out << '\n';
			}
			file_out << p(i,j) << ' ';
		}
	}
	file_out.close();
}

void NavierStokes2::SolveLaplacianP()
{
	int N(Ny+1), M(Nx+1);
	Matrix div(N,M);
	Div2thOrder(v, dx, dy, N, M, div);
	for(int i(0); i < N; ++i) {
		for(int j(0); j < M; ++j) {
			sec_membre_p(i,j) = (rho(i,j)*div(i,j))/dt + cond_bord_p(i,j);
		}
	}
	GradConjLaplacian2(dx, dy, N, M, 1e-12, N*M, p, sec_membre_p); // pression sous forme de vecteur
}

void NavierStokes2::Advance(int n, double tn)
{
	int N(Ny+1), M(Nx+1);
	precision a(dx*dx), b(dy*dy), y;
	// Actualisation des conditions aux bords (pour les conditions symétriques)
	
	Matrix gradpx(N,M), gradpy(N,M);
	Matrix pgz(N,M);

	const Matrix pconst(p),vgetX(v.GetAllVX()),vgetY(v.GetAllVY());



	if(firstIter) {
		firstIter = false;
	} else {
		spacescheme.SetVelocity(v);
		for(int i(0); i < N; ++i) {
			for(int j(0); j < M; ++j) {
				cond_bord_p(i,j) = 0.0; // initialisation
			}
			// Condition bords gauche et droite
			cond_bord_p(i,0) = -p.left(i,-1)/a;
			cond_bord_p(i,Nx) = -p.right(i,M)/a;
		}
		for(int j(0); j < M; ++j) {
			// Condition bords haut et bas 8
			cond_bord_p(0,j) += -p.bottom(-1,j)/b;
			cond_bord_p(Ny,j) += -p.top(N,j)/b;
		}
		spacescheme.SetD(D_nacl);
		timescheme_n.Advance(n, tn, spacescheme, zero);
		rho.Set(timescheme_n.GetIterate());
		SolveLaplacianP();
	}
	for(int i(0); i < N; ++i) {
		y = H-i*dy;
		for(int j(0); j < M; ++j) {
			pgz(i,j) = p(i,j)-rho(i,j)*g*y; // gravité
		}
	}
	pgz.CopyBoundaryCondition(p);
	Gradx2thOrder(dx, N, M, p, gradpx);
	Grady2thOrder(dy, N, M, pgz, gradpy);
	for(int i(0); i < N; ++i) {
		for(int j(0); j < M; ++j) {
			sec_membre_vx(i,j) = -(dt*gradpx(i,j))/rho(i,j);
			sec_membre_vy(i,j) = -(dt*gradpy(i,j))/rho(i,j);
		}
	}
	spacescheme.SetD(nu);
	timescheme_x.Advance(n, tn, spacescheme, sec_membre_vx);
	timescheme_y.Advance(n, tn, spacescheme, sec_membre_vy);
	v.SetAllVX(timescheme_x.GetIterate());
	v.SetAllVY(timescheme_y.GetIterate());
}




void GradConjLaplacian2(precision dx, precision dy, int N, int M, precision epsilon, int Nmax, Matrix& x, Matrix& b)
{
	int iter;
	Matrix d(N,M),w(N,M),r(N,M);
	precision alpha, beta, nr;

	Laplacian2thOrder(dx, dy, N, M, x, r);
	for(int i(0); i < N; ++i) {
		for(int j(0); j < M; ++j) {
			r(i,j) -= b(i,j);
			d(i,j) = r(i,j);
		}
	}
	iter = 1;
	nr = VecNorme(r);

	while(nr > epsilon && iter < Nmax) {
		Laplacian2thOrder(dx, dy, N, M, d, w);

		alpha = DotProduct(d,r)/DotProduct(d,w);

		for(int i(0); i < N; ++i) {
			for(int j(0); j < M; ++j) {
				x(i,j) -= alpha*d(i,j);
			}
		}

		beta = 1.0/(nr*nr);
		for(int i(0); i < N; ++i) {
			for(int j(0); j < M; ++j) {
				r(i,j) -= alpha*w(i,j);
			}
		}
		nr = VecNorme(r);
		beta = beta*nr*nr;
		for(int i(0); i < N; ++i) {
			for(int j(0); j < M; ++j) {
				d(i,j) = r(i,j)+beta*d(i,j);
			}
		}
		++iter;
	}
	cout << nr << ";\tnb_iter = " << iter << endl;
}

/**
 * Calcul du produit matrice vecteur : out = A.v
 * avec A la matrice du laplacien discrétisé à l'ordre 2
 */
void Laplacian2thOrder(precision dx, precision dy, int N, int M, const Matrix& v, Matrix& out)
{
	precision a(dx*dx), b(dy*dy);

	// 1er bloc
	out(0,0) = (v(0,1)-2.0*v(0,0))/a
					+ (v(1,0)-2.0*v(0,0))/b;
	for(int j(1); j < M-1; ++j) {
		out(0,j) = (v(0,j+1)-2.0*v(0,j)+v(0,j-1))/a
						+ (v(1,j)-2.0*v(0,j))/b;
	}
	out(0,M-1) = (v(Bij(0,M-2,M))-2.0*v(Bij(0,M-1,M)))/a
					  + (v(Bij(1,M-1,M))-2.0*v(Bij(0,M-1,M)))/b;

	// Milieu Matrice
	for(int i(1); i < N-1; ++i) {
		out(i,0) = (v(i,1)-2.0*v(i,0))/a
						+ (v(i+1,0)-2.0*v(i,0)+v(i-1,0))/b;
		for(int j(1); j < M-1; ++j) {
			out(i,j) = (v(i,j+1)-2.0*v(i,j)+v(i,j-1))/a
							+ (v(i+1,j)-2.0*v(i,j)+v(i-1,j))/b;
		}
		out(i,M-1) = (v(i,M-2)-2.0*v(i,M-1))/a
						  + (v(i+1,M-1)-2.0*v(i,M-1)+v(i-1,M-1))/b;
	}

	// Dernier bloc
	out(N-1,0) = (v(N-1,1)-2.0*v(N-1,0))/a
					+ (v(N-2,0)-2.0*v(N-1,0))/b;
	for(int j(1); j < M-1; ++j) {
		out(N-1,j) = (v(N-1,j+1)-2.0*v(N-1,j)+v(N-1,j-1))/a
						+ (v(N-2,j)-2.0*v(N-1,j))/b;
	}
	out(N-1,M-1) = (v(N-1,M-2)-2.0*v(N-1,M-1))/a
			+ (v(N-2,M-1)-2.0*v(N-1,M-1))/b;

}

void Div2thOrder(const Velocity& v, precision dx, precision dy, int N, int M, Matrix& out)
{
	precision a(2.0*dx), b(2.0*dy);


	// Milieu Matrice
	for(int i(0); i < N; ++i) {
		for(int j(0); j < M; ++j) {
			out(i,j) = (v.GetVX(i,j+1)-v.GetVX(i,j-1))/a
							+ (v.GetVY(i+1,j)-v.GetVY(i-1,j))/b;
		}
	}

}

void Gradx2thOrder(precision dx, int N, int M, const Matrix& v, Matrix& out)
{
	precision a(2.0*dx);

	for(int i(0); i < N; ++i) {
		for(int j(0); j < M; ++j) {
			out(i,j) = (v(i,j+1)-v(i,j-1))/a;
		}
	}
}



void Grady2thOrder(precision dy, int N, int M, const Matrix& v, Matrix& out)
{
	precision b(2.0*dy);

	for(int i(0); i < N; ++i) {
		for(int j(0); j < M; ++j) {
			out(i,j) = (v(i+1,j)-v(i-1,j))/b;
		}
	}
}

