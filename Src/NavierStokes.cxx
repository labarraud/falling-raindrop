/*
 * NavierStokes.cxx
 *
 *  Created on: 4 avr. 2016
 *      Author: tlariviere
 */

#include "../Include/NavierStokes.hxx"

NavierStokes::NavierStokes()
{ }

NavierStokes::NavierStokes(int _Nx,int _Ny,int _Nt,precision _L,precision _H,precision tfinal,Velocity& _v,Density& _rho,const Matrix& _p)
{
	SetInitialCondition(_Nx,_Ny,_Nt,_L,_H,tfinal,_v,_rho);
	_p.Mat2Vec(p);
}

void NavierStokes::SetInitialCondition(int _Nx,int _Ny,int _Nt,precision _L,precision _H,precision tfinal,Velocity& _v,Density& _rho)
{
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
	//nu = 1.007e-6; // viscosité cinématique de l'eau à 20°C
	nu = 5.0e-3;
	g = 9.81;
	int N(Ny+1), M(Nx+1), l(N*M);
	p.resize(l);
	cond_bord_p.resize(l);
	cond_bord_div_v.resize(l);
	cond_bord_gradx_p.resize(l);
	cond_bord_grady_p.resize(l);
	sec_membre_p.resize(l);
	sec_membre_vx.resize(l);
	sec_membre_vy.resize(l);
	timescheme_x.SetInitialCondition(0,dt,v.GetAllVY(),*this);
	timescheme_y.SetInitialCondition(0,dt,v.GetAllVX(),*this);
}

void NavierStokes::SetPressure(const Matrix& _p)
{
	_p.Mat2Vec(p);
}


void NavierStokes::WriteVtk(const string& nom) const
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
	file_out << "SCALARS velocityX float\n";
	file_out << "LOOKUP_TABLE default";
	for (int i = 0; i < N; i++) {
		for (int j=0; j<M;j++) {
			if((j%10)==0) {
				file_out << '\n';
			}
			file_out << v.GetVX(i,j) << ' ';
		}
	}
	file_out << "SCALARS velocityY float\n";
	file_out << "LOOKUP_TABLE default";
	for (int i = 0; i < N; i++) {
		for (int j=0; j<M;j++) {
			if((j%10)==0) {
				file_out << '\n';
			}
			file_out << v.GetVY(i,j) << ' ';
		}
	}
	file_out << "SCALARS pressure float\n";
	file_out << "LOOKUP_TABLE default";
	for (int i = 0; i < N; i++) {
		for (int j=0; j<M;j++) {
			if((j%10)==0) {
				file_out << '\n';
			}
			file_out << p[Bij(i,j,M)] << ' ';
		}
	}
	file_out.close();
}

void NavierStokes::SolveLaplacianP()
{
	int N(Ny+1), M(Nx+1), l(N*M);
	vector<precision> div(l);
	Div4thOrder(v, dx, dy, N, M, div);
	for(int i(0); i < N; ++i) {
		for(int j(0); j < M; ++j) {
			sec_membre_p[Bij(i,j,M)] = rho(i,j)*(div[Bij(i,j,M)]+cond_bord_div_v[Bij(i,j,M)])/dt + cond_bord_p[Bij(i,j,M)];
		}
	}
	GradConjLaplacian(dx, dy, N, M, 1e-6, N*M+1, p, sec_membre_p); // pression sous forme de vecteur
}

precision NavierStokes::UpwindY(double dt, double b, int i, int j, const Matrix& u)
{
	precision theta(dt/dy);
	int sign_b((b < 0) ? -1 : 1);
	precision u_ij(uij(i,j,u));
	return (u_ij-b*theta*sign_b*(3.0*uij(i+sign_b,j,u)+10.0*u_ij-18.0*uij(i-sign_b,j,u)+6.0*uij(i-2*sign_b,j,u)-uij(i-3*sign_b,j,u))/12.0);
}

precision NavierStokes::SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u)
{
	double sigma(dt/dx), ut;
	int sign_a((a < 0) ? -1 : 1);
	ut = UpwindY(dt, b, i, j, u);
	return (ut-a*sigma*sign_a*(3.0*UpwindY(dt,b,i,j+sign_a,u)+10.0*ut-18.0*UpwindY(dt,b,i,j-sign_a,u)+6.0*UpwindY(dt,b,i,j-2*sign_a,u)-UpwindY(dt,b,i,j-3*sign_a,u))/12.0);
}

void NavierStokes::AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y)
{

	//alpha = delta_t
	int N(Ny+1), M(Nx+1);
    for (int i=0; i<N; i++) {
    	for (int j=0; j<M; j++) {
			y(i,j) += SplittingX(alpha, v.GetVX(i,j), v.GetVY(i,j), i, j, u) - u(i,j)
						+ alpha*nu*(16.0*uij(i+1,j, u)-uij(i+2,j, u)-30.0*u(i,j)-uij(i-2,j, u)+16.0*uij(i-1,j, u))/(dx*dx)
						+ alpha*nu*(16.0*uij(i,j+1, u)-uij(i,j+2, u)-30.0*u(i,j)-uij(i,j-2, u)+16.0*uij(i,j-1, u))/(dy*dy);
		}
   }
}

void NavierStokes::Advance(int n, double tn)
{
	int N(Ny+1), M(Nx+1);
	precision a1(12.0*dx*dx), b1(12.0*dy*dy), a2(12.0*dx), b2(12.0*dy);
	// Actualisation des conditions aux bords (pour les conditions symétriques)
	for(int i(0); i < N; ++i) {
		for(int j(0); j < M; ++j) {
			cond_bord_p[Bij(i,j,M)] = 0.0; // initialisation
			cond_bord_div_v[Bij(i,j,M)] = 0.0;
			cond_bord_gradx_p[Bij(i,j,M)] = 0.0;
			cond_bord_grady_p[Bij(i,j,M)] = 0.0;
		}
		// Condition bords gauche et droite
		cond_bord_p[Bij(i,0,M)] = -15.0*p_bord_gauche(i,-2)/a1;
		cond_bord_p[Bij(i,1,M)] = p_bord_gauche(i,-1)/a1;
		cond_bord_p[Bij(i,Nx-1,M)] = p_bord_droit(i,M)/a1;
		cond_bord_p[Bij(i,Nx,M)] = -15.0*p_bord_droit(i,M+1)/a1;

		cond_bord_gradx_p[Bij(i,0,M)] = 7.0*p_bord_gauche(i,-2)/a2;
		cond_bord_gradx_p[Bij(i,1,M)] = -p_bord_gauche(i,-1)/a2;
		cond_bord_gradx_p[Bij(i,Nx-1,M)] = p_bord_droit(i,M)/a2;
		cond_bord_gradx_p[Bij(i,Nx,M)] = -7.0*p_bord_droit(i,M+1)/a2;

		cond_bord_div_v[Bij(i,0,M)] = 7.0*vx_bord_gauche(i,-2,v.GetAllVX())/a2;
		cond_bord_div_v[Bij(i,1,M)] = -vx_bord_gauche(i,-1,v.GetAllVX())/a2;
		cond_bord_div_v[Bij(i,Nx-1,M)] = vx_bord_droit(i,M,v.GetAllVX())/a2;
		cond_bord_div_v[Bij(i,Nx,M)] = -7.0*vx_bord_droit(i,M+1,v.GetAllVX())/a2;
	}
	for(int j(0); j < M; ++j) {
		// Condition bords haut et bas
		cond_bord_p[Bij(0,j,M)] += -15.0*p_bord_haut(-2,j)/b1;
		cond_bord_p[Bij(1,j,M)] += p_bord_haut(-1,j)/b1;
		cond_bord_p[Bij(Ny-1,j,M)] += p_bord_bas(N,j)/b1;
		cond_bord_p[Bij(Ny,j,M)] += -15.0*p_bord_bas(N+1,j)/b1;

		cond_bord_grady_p[Bij(0,j,M)] = 7.0*p_bord_haut(-2,j)/b2;
		cond_bord_grady_p[Bij(1,j,M)] = -p_bord_haut(-1,j)/b2;
		cond_bord_grady_p[Bij(Ny-1,j,M)] = p_bord_bas(N,j)/b2;
		cond_bord_grady_p[Bij(Ny,j,M)] = -7.0*p_bord_bas(N+1,j)/b2;

		cond_bord_div_v[Bij(0,j,M)] += 7.0*vy_bord_haut(-2,j,v.GetAllVY())/b2;
		cond_bord_div_v[Bij(1,j,M)] += -vy_bord_haut(-1,j,v.GetAllVY())/b2;
		cond_bord_div_v[Bij(Ny-1,j,M)] += vy_bord_bas(N,j,v.GetAllVY())/b2;
		cond_bord_div_v[Bij(Ny,j,M)] += -7.0*vy_bord_bas(N+1,j,v.GetAllVY())/b2;
	}
	vector<precision> gradpx(N*M), gradpy(N*M);

	SolveLaplacianP();
	Gradx4thOrder(dx, N, M, p, gradpx);
	Grady4thOrder(dy, N, M, p, gradpy);
	for(int i(0); i < N; ++i) {
		for(int j(0); j < M; ++j) {
			sec_membre_vx[Bij(i,j,M)] = -(dt*(gradpx[Bij(i,j,M)]+cond_bord_gradx_p[Bij(i,j,M)]))/rho(i,j);
			sec_membre_vy[Bij(i,j,M)] = -(dt*(gradpy[Bij(i,j,M)]+cond_bord_grady_p[Bij(i,j,M)]))/rho(i,j);//-dt*g;
		}
	}
	vx = true;
	timescheme_x.Advance(n, tn, *this);
	vx = false;
	timescheme_y.Advance(n, tn, *this);
	v.SetAllVX(timescheme_x.GetIterate());
	v.SetAllVY(timescheme_y.GetIterate());
}

precision NavierStokes::uij(int i, int j, const Matrix& u)
{
	if(i < 0) {
		return (vx ? vx_bord_haut(i,j,u) : vy_bord_haut(i,j,u));
	} else if(i > Ny) {
		return (vx ? vx_bord_bas(i,j,u) : vy_bord_bas(i,j,u));
	} else if(j < 0) {
		return (vx ? vx_bord_gauche(i,j,u) : vy_bord_gauche(i,j,u));
	} else if(j > Nx) {
		return (vx ? vx_bord_droit(i,j,u) : vy_bord_droit(i,j,u));
	} else {
		return u(i,j);
	}
}

precision NavierStokes::p_bord_droit(int i, int j) const { return 0.0; }
precision NavierStokes::p_bord_gauche(int i, int j) const { return 0.0; }
precision NavierStokes::p_bord_haut(int i, int j) const { return 0.0; }
precision NavierStokes::p_bord_bas(int i, int j) const { return 0.0; }
precision NavierStokes::vx_bord_droit(int i, int j, const Matrix& u) const {
	int M(u.GetM());
	if(j == M) {
		return 8.0*(u(i,M-1)-u(i,M-3))+u(i,M-4);
	} else if( j == M+1) {
		return 8.0*u(i,M-2)-u(i,M-3);
	} else {
		return 0.0;
	}
}

precision NavierStokes::vx_bord_gauche(int i, int j, const Matrix& u) const { return 1.0; }
precision NavierStokes::vx_bord_haut(int i, int j, const Matrix& u) const { return 0.0; }
precision NavierStokes::vx_bord_bas(int i, int j, const Matrix& u) const { return 0.0; }
precision NavierStokes::vy_bord_droit(int i, int j, const Matrix& u) const { return 0.0; }
precision NavierStokes::vy_bord_gauche(int i, int j, const Matrix& u) const { return 0.0; }
precision NavierStokes::vy_bord_haut(int i, int j, const Matrix& u) const { return 0.0; }
precision NavierStokes::vy_bord_bas(int i, int j, const Matrix& u) const { return 0.0; }



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
	precision a(12.0*dx), b(12.0*dy);


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
					+ (v.GetVY(N-3,0)-8.0*v.GetVY(N-2,0))/b;
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

void Gradx4thOrder(precision dx, int N, int M, const vector<precision>& v, vector<precision>& out)
{
	precision a(12.0*dx);

	for(int i(0); i < N; ++i) {
		out[Bij(i,0,M)] = (8.0*v[Bij(i,1,M)]-v[Bij(i,2,M)])/a;
		out[Bij(i,1,M)] = (8.0*v[Bij(i,2,M)]-v[Bij(i,3,M)]-8.0*v[Bij(i,0,M)])/a;
		for(int j(2); j < M-2; ++j) {
			out[Bij(i,j,M)] = (8.0*v[Bij(i,j+1,M)]-v[Bij(i,j+2,M)]+v[Bij(i,j-2,M)]-8.0*v[Bij(i,j-1,M)])/a;
		}
		out[Bij(i,M-2,M)] = (8.0*v[Bij(i,M-1,M)]+v[Bij(i,M-4,M)]-8.0*v[Bij(i,M-3,M)])/a;
		out[Bij(i,M-1,M)] = (v[Bij(i,M-3,M)]-8.0*v[Bij(i,M-2,M)])/a;
	}
}



void Grady4thOrder(precision dy, int N, int M, const vector<precision>& v, vector<precision>& out)
{
	precision b(12.0*dy);

	for(int j(0); j < M; ++j) {
		out[Bij(0,j,M)] = (8.0*v[Bij(1,j,M)]-v[Bij(2,j,M)])/b;
		out[Bij(1,j,M)] = (8.0*v[Bij(2,j,M)]-v[Bij(3,j,M)]-8.0*v[Bij(0,j,M)])/b;
		for(int i(2); i < N-2; ++i) {
			out[Bij(i,j,M)] = (8.0*v[Bij(i+1,j,M)]-v[Bij(i+2,j,M)]+v[Bij(i-2,j,M)]-8.0*v[Bij(i-1,j,M)])/b;
		}
		out[Bij(N-2,j,M)] = (8.0*v[Bij(N-1,j,M)]+v[Bij(N-4,j,M)]-8.0*v[Bij(N-3,j,M)])/b;
		out[Bij(N-1,j,M)] = (v[Bij(N-3,j,M)]-8.0*v[Bij(N-2,j,M)])/b;
	}
}

