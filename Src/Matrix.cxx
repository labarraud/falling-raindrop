#ifndef MATRIX_CXX
#define MATRIX_CXX

#include "../Include/Matrix.hxx"
//test

Matrix::Matrix()
	:	N(0),M(0)
{ }

Matrix::Matrix(const Matrix& m)
	:	N(m.GetN()),M(m.GetM())
{
	val.resize(N);
	for(int i(0); i < N; ++i) {
		val[(unsigned)i].resize((unsigned)M);
		for(int j(0); j < M; ++j) {
			//val[(unsigned)i][(unsigned)j]
				(*this)(i,j) = m(i,j);
		}
	}
}

Matrix::Matrix(int _N, int _M)
:	N(_N),M(_M)
{
	val.resize(N);
	for(int i(0); i < N; ++i) {
		val[(unsigned)i].resize((unsigned)M);
		for(int j(0); j < M; ++j) {
				(*this)(i,j) = 0;
		}
	}
}

Matrix::Matrix(int _N)
	:	N(_N),M(_N)
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
	return val[(unsigned)i][(unsigned)j];
}

const precision& Matrix::operator()(int i, int j) const
{
	return val[(unsigned)i][(unsigned)j];
}

Matrix& Matrix::operator=(const Matrix& m)
{
	N = m.GetN();
	M = m.GetM();
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
	Matrix var(m.GetM(),m.GetN());
	int N = m.GetN();
	int M = m.GetM();
	for(int i(0); i < N; ++i) {
		for (int j(0); j<M; ++j)
			var(i,j) = a*m(i,j);
		}
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
	return var/(N*M);
}

void Matrix::GradConj(const Matrix& A, Matrix P, Matrix X)
{


#endif
