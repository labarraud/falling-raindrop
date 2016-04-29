#ifndef MATRIX_CXX
#define MATRIX_CXX

#include <iostream>
#include <string>
#include <cmath>
#include "../Include/Matrix.hxx"
//test

Matrix::Matrix()
	:	N(0),M(0),X0(0),Xn(0),Y0(0),Yn(0),bcX0(periodic),bcXn(periodic),bcY0(periodic),bcYn(periodic)
{ }

Matrix::Matrix(const Matrix& m)
	:	N(m.GetN()),M(m.GetM())
{
	val.resize((unsigned)N);
	for(int i(0); i < N; ++i) {
		val[(unsigned)i].resize((unsigned)M);
		for(int j(0); j < M; ++j) {
			//val[(unsigned)i][(unsigned)j]
				(*this)(i,j) = m(i,j);
		}
	}

	this->X0=m.X0;
	this->Xn=m.Xn;
	this->Y0=m.Y0;
	this->Yn=m.Yn;

	this->bcX0=m.bcX0;
	this->bcXn=m.bcXn;
	this->bcY0=m.bcY0;
	this->bcYn=m.bcYn;
}

Matrix::Matrix(int _N, int _M)
:	N(_N),M(_M),X0(0),Xn(0),Y0(0),Yn(0),bcX0(periodic),bcXn(periodic),bcY0(periodic),bcYn(periodic)
{
	val.resize((unsigned)N);
	for(int i(0); i < N; ++i) {
		val[(unsigned)i].resize((unsigned)M);
		for(int j(0); j < M; ++j) {
				(*this)(i,j) = 0;
		}
	}
}

Matrix::Matrix(int _N)
	:	N(_N),M(_N),X0(0),Xn(0),Y0(0),Yn(0),bcX0(periodic),bcXn(periodic),bcY0(periodic),bcYn(periodic)
{
	val.resize(N);
		for(int i(0); i < N; ++i) {
			val[(unsigned)i].resize((unsigned)N);
			for(int j(0); j < N; ++j) {
					(*this)(i,j) = 0;
			}
	}

}

int Matrix::GetN() const
{
	return N;
}

int Matrix::GetM() const
{
	return M;
}	

void Matrix::Zero()
{
	for(int i(0); i < N; ++i) {
		for(int j(0); j < M; ++j) {
				(*this)(i,j) = 0;
		}
	}
}

void Matrix::SetBoundaryCondition(BoundaryCondition bcX0,precision X0,BoundaryCondition bcXn,precision Xn, BoundaryCondition bcY0, precision Y0,BoundaryCondition bcYn, precision Yn)
{
	this->X0=X0;
	this->Xn=Xn;
	this->Y0=Y0;
	this->Yn=Yn;

	this->bcX0=bcX0;
	this->bcXn=bcXn;
	this->bcY0=bcY0;
	this->bcYn=bcYn;

}

BoundaryCondition const Matrix::getbcX0() const { return bcX0; }
BoundaryCondition const Matrix::getbcXn() const { return bcXn; }
BoundaryCondition const Matrix::getbcY0() const { return bcY0; }
BoundaryCondition const Matrix::getbcYn() const { return bcYn; }
precision const Matrix::getX0() const { return X0; }
precision const Matrix::getXn() const { return Xn; }
precision const Matrix::getY0() const { return Y0; }
precision const Matrix::getYn() const { return Yn; }



void Matrix::Reallocate(int _N)
{
	N = _N;
	M = _N;
	val.resize(N);
	for(int i(0); i < N; ++i) {
		val[(unsigned)i].resize((unsigned)N);
	}
}

void Matrix::Reallocate(int _N, int _M)
{
	N = _N;
	M = _M;
	val.resize((unsigned)N);
	for(int i(0); i < N; ++i) {
		val[(unsigned)i].resize((unsigned)M);
	}
}

void Matrix::Clear()
{


	for(int i(0); i < N; ++i) {
		val[(unsigned)i].clear();
	}
	val.clear();
	N = 0;
	M = 0;
}



precision& Matrix::operator()(int i, int j)
{

	if(i>N || j>M) {
		cout << "i=" << i << "; N=" << N << "; j=" << j << "; M=" << M << endl;
	}
	return val[(unsigned)i][(unsigned)j];

}



const precision Matrix::bottom(int i, int j) const
{
	int ivar=i, jvar=j;
	switch (this->bcY0)
	{
		case periodic:
			ivar=N+i;
			if (j<0)
				jvar=M+j;
			else if (j>N)
				jvar=j%M;
			return val[(unsigned)ivar][(unsigned)jvar];
			break;
		case dirichlet:
			return this->Y0;
			break;
		case neumann:
			cout << "neumann bottom not define"<< endl;
			return 0;
			break;
	}
	return 0;
}

const precision Matrix::top(int i, int j) const
{
	int ivar=i, jvar=j;

	switch (this->bcYn)
	{
		case periodic:
			ivar=i%N;
			if (j<0)
				jvar=M+j;
			else if (j>N)
				jvar=j%M;
			return val[(unsigned)ivar][(unsigned)jvar];
			break;
		case dirichlet:
			return this->Yn;
			break;
		case neumann:
			cout << "neumann top not define"<< endl;
			return 0;
			break;
	}
	return 0;
}


const precision Matrix::left(int i, int j) const
{
	int ivar, jvar;

	switch (this->bcX0)
	{
		case periodic:
			ivar=i;
			jvar=M+j;
			return val[(unsigned)ivar][(unsigned)jvar];
			break;
		case dirichlet:
			return this->X0;
			break;
		case neumann:
			cout << "neumann left not define"<< endl;
			return 0;
			break;
	}
	return 0;
}

const precision Matrix::right(int i, int j) const
{
	int ivar, jvar;

	switch (this->bcXn)
	{
		case periodic:
			ivar=i;
			jvar=j%M;
			return val[(unsigned)ivar][(unsigned)jvar];
			break;
		case dirichlet:
			return this->Xn;
			break;
		case neumann:
			cout << "neumann right not define"<< endl;
			return 0;
			break;
	}
	return 0;
}


const precision Matrix::operator()(int i, int j) const
{

		if(i < 0) {
			return bottom(i,j);
		} else if(i >= N) {
			return top(i,j);
		} else if(j < 0) {
			return left(i,j);
		} else if(j >= M) {
			return right(i,j);
		} else {
			return val[(unsigned)i][(unsigned)j];
		}



}

Matrix& Matrix::operator=(const Matrix& m)
{
	N = m.GetN();
	M = m.GetM();

	this->X0=m.X0;
	this->Xn=m.Xn;
	this->Y0=m.Y0;
	this->Yn=m.Yn;

	this->bcX0=m.bcX0;
	this->bcXn=m.bcXn;
	this->bcY0=m.bcY0;
	this->bcYn=m.bcYn;
	
	this->val.resize(N);
	for(int i(0); i < N; ++i) {
	  val[(unsigned)i].resize(M);
	  for (int j(0); j<M; ++j)
	    (*this)(i,j) = m(i,j);
	}
	return *this;
}

Matrix Matrix::operator+(const Matrix& m)
{
	Matrix var(m.GetN(),m.GetM());
	for(int i(0); i < N; ++i) {
		for (int j(0); j<M; ++j)
			var(i,j) = (*this)(i,j) + m(i,j);
	}
	var.X0=this->X0;
	var.Xn=this->Xn;
	var.Y0=this->Y0;
	var.Yn=this->Yn;

	var.bcX0=this->bcX0;
	var.bcXn=this->bcXn;
	var.bcY0=this->bcY0;
	var.bcYn=this->bcYn;

	return var;
}

ostream& operator<<(ostream& out, const Matrix& m)
{
	int N(m.GetN());
	int M(m.GetM());
	for(int i(0); i < N; ++i) {
		for(int j = 0; j < M; ++j) {
			out << m(i,j) << "\t";
		}
		out << endl;
	}
	return out;
}

Matrix operator*(precision a, const Matrix& m)
{
	int N = m.GetN();
	int M = m.GetM();
	Matrix var(N,M);
	for(int i(0); i < N; ++i) {
		for (int j(0); j<M; ++j)
			var(i,j) = a*m(i,j);
		}

	var.SetBoundaryCondition(m.getbcX0(),m.getX0(),m.getbcXn(),m.getXn(), m.getbcY0(), m.getY0(),m.getbcYn(), m.getYn());
	return var;


}

precision Matrix::norme2()
{
	precision var=0;
	for(int i(0); i < N; ++i)
	{
			for (int j(0); j<M; ++j)
				var=var+(*this)(i,j)*(*this)(i,j);
	}
	return sqrt(var);
}

precision Matrix::distnorme2(const Matrix& m)
{
	precision var=0;
	for(int i(0); i < N; ++i)
	{
			for (int j(0); j<M; ++j)
				var=var +((*this)(i,j)-m(i,j))*((*this)(i,j)-m(i,j));
	}
	return sqrt(var);
}

void Matrix::WriteGnuPlot(const string& nom) const
	{
		ofstream file_out(nom.data());
		file_out.precision(15);
	//	double x,y;
		for (int i = 0; i < N; i++)
			{
			for (int j=0; j<M;j++)
			{
				file_out << (*this)(i,j)<< " ";
			}
				file_out << "\n";
		  }
		  file_out.close();
	}




void Matrix::WriteVtk(const string& nom, precision dx, precision dy) const
	{
		  ofstream file_out(nom.data());
		  file_out << "# vtk DataFile Version 2.0\n";
		  file_out << "Titre\n";
		  file_out << "ASCII\n";
		  file_out << "DATASET STRUCTURED_POINTS\n";
		  file_out << "DIMENSIONS " << M << " " << N << " 1\n";
		  file_out << "ORIGIN 0.0 0.0 0.0\n";
		  file_out << "SPACING " << dx << " " << dy << " 0.0\n";
		  file_out << "POINT_DATA " << (M*N) << " \n";
		  file_out << "SCALARS rho float\n";
		  file_out << "LOOKUP_TABLE default";
		  for (int i = 0; i < N; i++) {
			for (int j=0; j<M;j++) {
			  	if((j%10)==0) {
			  		file_out << '\n';
			  	}
				file_out << (*this)(i,j) << ' ';
			}
		  }
		  file_out.close();
	}

void Matrix::Mat2Vec(vector<precision>& out) const
{
	out.resize(N*M);
	for(int i(0),j; i < N; ++i) {
		for(j = 0; j < M; ++j) {
			out[unsigned(i*M+j)]=(*this)(i,j);
		}
	}
}

void Matrix::Vec2Mat(const vector<precision>& in)
{
	for(int i(0),j; i < N; ++i) {
		for(j = 0; j < M; ++j) {
			(*this)(i,j)=in[unsigned(i*M+j)];
		}
	}
}


precision VecNorme(const vector<precision>& v)
{
	precision S(0.0);
	for(unsigned int i(0),l(v.size()); i < l; ++i) {
		S += v[i]*v[i];
	}
	return sqrt(S);
}

precision DotProduct(const vector<precision>& v1, const vector<precision>& v2)
{
	precision S(0.0);
	for(unsigned int i(0),l(v1.size()); i < l; ++i) {
		S += v1[i]*v2[i];
	}
	return S;
}

/**
 * Bijection coordonnées matrice de taille N*M : (i,j) vers coordonnées vecteur
 * (stockage ligne par ligne de la matrice)
 */
unsigned int Bij(int i, int j, int M)
{
	return unsigned(M*i+j);
}

#endif
