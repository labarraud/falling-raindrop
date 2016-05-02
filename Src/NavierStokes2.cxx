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
	rho_mer = 1032.0;
	p_atm = 1015000.0;
	g = 9.81;
	int N(Ny+1), M(Nx+1), l(N*M);
	cond_bord_p.resize(l);
	cond_bord_div_v.resize(l);
	cond_bord_gradx_p.resize(l);
	cond_bord_grady_p.resize(l);
	sec_membre_p.resize(l);
	sec_membre_vx.resize(l);
	sec_membre_vy.resize(l);
	timescheme_x.SetInitialCondition(0,dt,v.GetAllVX(),*this);
	timescheme_y.SetInitialCondition(0,dt,v.GetAllVY(),*this);
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
			file_out << p(Bij(i,j,M)) << ' ';
		}
	}
	file_out.close();
}

void NavierStokes2::SolveLaplacianP()
{
	int N(Ny+1), M(Nx+1), l(N*M);
	vector<precision> div(l);
	Div2thOrder(v, dx, dy, N, M, div);
	for(int i(0); i < N; ++i) {
		for(int j(0); j < M; ++j) {
			sec_membre_p[Bij(i,j,M)] = rho(i,j)*(div[Bij(i,j,M)]+cond_bord_div_v[Bij(i,j,M)])/dt + cond_bord_p[Bij(i,j,M)];
		}
	}
	GradConjLaplacian2(dx, dy, N, M, 1e-12, l, p, sec_membre_p); // pression sous forme de vecteur
}

precision NavierStokes2::UpwindY(double dt, double b, int i, int j, const Matrix& u)
{
	precision theta(dt/dy);
	int sign_b((b < 0) ? -1 : 1);
	precision u_ij(u(i,j));
	return (u_ij - b*theta*sign_b*(3.0*u_ij-4.0*u(i-sign_b,j)+u(i-2*sign_b,j))/2.0);
}

precision NavierStokes2::SplittingX(precision dt, precision a, precision b, int i, int j, const Matrix& u)
{
	double sigma(dt/dx), ut;
	int sign_a((a < 0) ? -1 : 1);
	ut = UpwindY(dt, b, i, j, u);
	return (ut - a*sigma*sign_a*(3.0*ut-4.0*UpwindY(dt,b,i,j-sign_a,u)+UpwindY(dt,b,i,j-2*sign_a,u))/2.0);
}

void NavierStokes2::AddFunction(precision alpha, const Matrix& u, precision t, Matrix& y, const vector<precision>& sec_membre)
{

	//alpha = delta_t
	int N(Ny+1), M(Nx+1);
	precision a(dx*dx), b(dy*dy);
    for (int i=0; i<N; i++) {
    	for (int j=0; j<M; j++) {
			y(i,j) += SplittingX(alpha, v.GetVX(i,j), v.GetVY(i,j), i, j, u) - u(i,j)
						+ alpha*nu*(u(i,j+1)-2.0*u(i,j)+u(i,j-1))/a
						+ alpha*nu*(u(i+1,j)-2.0*u(i,j)+u(i-1,j))/b
						+ sec_membre[Bij(i,j,M)];
		}
   }
}

void NavierStokes2::Advance(int n, double tn)
{
	int N(Ny+1), M(Nx+1);
	precision a1(dx*dx), b1(dy*dy), a2(2.0*dx), b2(2.0*dy), y;
	// Actualisation des conditions aux bords (pour les conditions symétriques)
	
	vector<precision> gradpx(N*M), gradpy(N*M);
	Matrix pgz(N,M);

	const Matrix pconst(p),vgetX(v.GetAllVX()),vgetY(v.GetAllVY());



	if(firstIter) {
		firstIter = false;
	} else {
		for(int i(0); i < N; ++i) {
			for(int j(0); j < M; ++j) {
				cond_bord_p[Bij(i,j,M)] = 0.0; // initialisation
				cond_bord_div_v[Bij(i,j,M)] = 0.0;
			}
			// Condition bords gauche et droite
			cond_bord_p[Bij(i,0,M)] = -pconst(i,-1)/a1;
			cond_bord_p[Bij(i,Nx,M)] = -pconst(i,M)/a1;

			cond_bord_div_v[Bij(i,0,M)] = -vgetY(i,-1)/a2;
			cond_bord_div_v[Bij(i,Nx,M)] = vgetY(i,M)/a2;
		}
		for(int j(0); j < M; ++j) {
			// Condition bords haut et bas
			cond_bord_p[Bij(0,j,M)] += -pconst(-1,j)/b1;
			cond_bord_p[Bij(Ny,j,M)] += -pconst(N,j)/b1;

			cond_bord_div_v[Bij(0,j,M)] += -vgetY(-1,j)/b2;
			cond_bord_div_v[Bij(Ny,j,M)] += vgetY(N,j)/b2;
		}
		SolveLaplacianP();
	}
	for(int i(0); i < N; ++i) {
		y = H-i*dy;
		for(int j(0); j < M; ++j) {
			cond_bord_gradx_p[Bij(i,j,M)] = 0.0;
			cond_bord_grady_p[Bij(i,j,M)] = 0.0;
			pgz(Bij(i,j,M)) = p(Bij(i,j,M))-rho(i,j)*g*y; // gravité
		}
		// Condition bords gauche et droite

		cond_bord_gradx_p[Bij(i,0,M)] = -pconst(i,-1)/a2;
		cond_bord_gradx_p[Bij(i,Nx,M)] = pconst(i,M)/a2;
	}
	for(int j(0); j < M; ++j) {
		// Condition bords haut et bas
		cond_bord_grady_p[Bij(0,j,M)] = -pconst(-1,j)/b2;
		cond_bord_grady_p[Bij(Ny,j,M)] = pconst(N,j)/b2;
	}
	Gradx2thOrder(dx, N, M, p, gradpx);
	Grady2thOrder(dy, N, M, pgz, gradpy);
	for(int i(0); i < N; ++i) {
		for(int j(0); j < M; ++j) {
			sec_membre_vx[Bij(i,j,M)] = -(dt*(gradpx[Bij(i,j,M)]+cond_bord_gradx_p[Bij(i,j,M)]))/rho(i,j);
			sec_membre_vy[Bij(i,j,M)] = -(dt*(gradpy[Bij(i,j,M)]+cond_bord_grady_p[Bij(i,j,M)]))/rho(i,j);
		}
	}
	timescheme_x.Advance(n, tn, *this, sec_membre_vx);
	timescheme_y.Advance(n, tn, *this, sec_membre_vy);
	v.SetAllVX(timescheme_x.GetIterate());
	v.SetAllVY(timescheme_y.GetIterate());
}




void GradConjLaplacian2(precision dx, precision dy, int N, int M, precision epsilon, int Nmax, Matrix& x, const vector<precision>& b)
{
	unsigned int l(N*M);
	int iter;
	vector<precision> w(l),r(l);
	Matrix d(N,M);
	precision alpha, beta, nr;

	Laplacian2thOrder(dx, dy, N, M, x,r);
	for(int i(0); i < l; ++i) {
		r[(unsigned)i] -= b[i];
		d(i)=r[(unsigned)i];
	}
	iter = 1;
	nr = VecNorme(r);

	while(nr > epsilon && iter < Nmax) {
		Laplacian2thOrder(dx, dy, N, M, d,w);

		alpha = DotProduct(d,r)/DotProduct(d,w);

		for(int i(0); i < l; ++i) {
			x(i) -= alpha*d(i);
		}

		beta = 1.0/(nr*nr);
		for(int i(0); i < l; ++i) {
			r[(unsigned)i] -= alpha*w[(unsigned)i];
		}
		nr = VecNorme(r);
		beta = beta*nr*nr;
		for( int i(0); i < l; ++i) {
			d(i) = r[(unsigned)i]+beta*d(i);
		}
		cout << nr << endl;
		++iter;
	}
	cout << nr << ";\tnb_iter = " << iter << endl;
}

/**
 * Calcul du produit matrice vecteur : out = A.v
 * avec A la matrice du laplacien discrétisé à l'ordre 2
 */
void Laplacian2thOrder(precision dx, precision dy, int N, int M, const Matrix& v, vector<precision>& out)
{
	precision a(dx*dx), b(dy*dy);

	// 1er bloc
	out[Bij(0,0,M)] = (v(Bij(0,1,M))-2.0*v(Bij(0,0,M)))/a
					+ (v(Bij(1,0,M))-2.0*v(Bij(0,0,M)))/b;
	for(int j(1); j < M-1; ++j) {
		out[Bij(0,j,M)] = (v(Bij(0,j+1,M))-2.0*v(Bij(0,j,M))+v(Bij(0,j-1,M)))/a
						+ (v(Bij(1,j,M))-2.0*v(Bij(0,j,M)))/b;
	}
	out[Bij(0,M-1,M)] = (v(Bij(0,M-2,M))-2.0*v(Bij(0,M-1,M)))/a
					  + (v(Bij(1,M-1,M))-2.0*v(Bij(0,M-1,M)))/b;

	// Milieu Matrice
	for(int i(1); i < N-1; ++i) {
		out[Bij(i,0,M)] = (v(Bij(i,1,M))-2.0*v(Bij(i,0,M)))/a
						+ (v(Bij(i+1,0,M))-2.0*v(Bij(i,0,M))+v(Bij(i-1,0,M)))/b;
		for(int j(1); j < M-1; ++j) {
			out[Bij(i,j,M)] = (v(Bij(i,j+1,M))-2.0*v(Bij(i,j,M))+v(Bij(i,j-1,M)))/a
							+ (v(Bij(i+1,j,M))-2.0*v(Bij(i,j,M))+v(Bij(i-1,j,M)))/b;
		}
		out[Bij(i,M-1,M)] = (v(Bij(i,M-2,M))-2.0*v(Bij(i,M-1,M)))/a
						  + (v(Bij(i+1,M-1,M))-2.0*v(Bij(i,M-1,M))+v(Bij(i-1,M-1,M)))/b;
	}

	// Dernier bloc
	out[Bij(N-1,0,M)] = (v(Bij(N-1,1,M))-2.0*v(Bij(N-1,0,M)))/a
					+ (v(Bij(N-2,0,M))-2.0*v(Bij(N-1,0,M)))/b;
	for(int j(1); j < M-1; ++j) {
		out[Bij(N-1,j,M)] = (v(Bij(N-1,j+1,M))-2.0*v(Bij(N-1,j,M))+v(Bij(N-1,j-1,M)))/a
						+ (v(Bij(N-2,j,M))-2.0*v(Bij(N-1,j,M)))/b;
	}
	out[Bij(N-1,M-1,M)] = (v(Bij(N-1,M-2,M))-2.0*v(Bij(N-1,M-1,M)))/a
			+ (v(Bij(N-2,M-1,M))-2.0*v(Bij(N-1,M-1,M)))/b;

}

void Div2thOrder(const Velocity& v, precision dx, precision dy, int N, int M, vector<precision>& out)
{
	precision a(2.0*dx), b(2.0*dy);


	// 1er bloc
	out[Bij(0,0,M)] = v.GetVX(0,1)/a
					+ v.GetVY(1,0)/b;
	for(int j(1); j < M-1; ++j) {
		out[Bij(0,j,M)] = (v.GetVX(0,j+1)-v.GetVX(0,j-1))/a
						+ v.GetVY(1,j)/b;
	}
	out[Bij(0,M-1,M)] = -v.GetVX(0,M-2)/a
					+ v.GetVY(1,M-1)/b;

	// Milieu Matrice
	for(int i(1); i < N-1; ++i) {
		out[Bij(i,0,M)] = v.GetVX(i,1)/a
						+ (v.GetVY(i+1,0)-v.GetVY(i-1,0))/b;
		for(int j(1); j < M-1; ++j) {
			out[Bij(i,j,M)] = (v.GetVX(i,j+1)-v.GetVX(i,j-1))/a
							+ (v.GetVY(i+1,j)-v.GetVY(i-1,j))/b;
		}
		out[Bij(i,M-1,M)] = -v.GetVX(i,M-2)/a
						+ (v.GetVY(i+1,M-1)-v.GetVY(i-1,M-1))/b;
	}

	// Dernier bloc
	out[Bij(N-1,0,M)] = v.GetVX(N-1,1)/a
					- v.GetVY(N-2,0)/b;
	for(int j(1); j < M-1; ++j) {
		out[Bij(N-1,j,M)] = (v.GetVX(N-1,j+1)-v.GetVX(N-1,j-1))/a
						- v.GetVY(N-2,j)/b;
	}
	out[Bij(N-1,M-1,M)] = -v.GetVX(N-1,M-2)/a
					- v.GetVY(N-2,M-1)/b;

}

void Gradx2thOrder(precision dx, int N, int M, const Matrix& v, vector<precision>& out)
{
	precision a(2.0*dx);

	for(int i(0); i < N; ++i) {
		out[Bij(i,0,M)] = v(Bij(i,1,M))/a;
		for(int j(1); j < M-1; ++j) {
			out[Bij(i,j,M)] = (v(Bij(i,j+1,M))-v(Bij(i,j-1,M)))/a;
		}
		out[Bij(i,M-1,M)] = -v(Bij(i,M-2,M))/a;
	}
}



void Grady2thOrder(precision dy, int N, int M, const Matrix& v, vector<precision>& out)
{
	precision b(2.0*dy);

	for(int j(0); j < M; ++j) {
		out[Bij(0,j,M)] = v(Bij(1,j,M))/b;
		for(int i(1); i < N-1; ++i) {
			out[Bij(i,j,M)] = (v(Bij(i+1,j,M))-v(Bij(i-1,j,M)))/b;
		}
		out[Bij(N-1,j,M)] = -v(Bij(N-2,j,M))/b;
	}
}

