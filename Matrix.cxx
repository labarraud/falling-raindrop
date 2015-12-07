#ifndef MATRIX_CXX
#define MATRIX_CXX

#include "Matrix.hxx"


Matrix::Matrix()
	:	N(0)
{ }

Matrix::Matrix(const Matrix& m)
	:	N(m.GetM())
{
	for(int i(0); i < N; ++i) {
		val(i) = m.val(i);
	}
}

Matrix::Matrix(int _N)
	:	N(_N), val(N)
{
	for(int i(0); i < N; ++i) {
		val(i).Reallocate(N);
	}
}

int Matrix::GetM() const
{
	return N;
}	

void Matrix::Zero()
{
	for(int i(0); i < N; ++i) {
		val(i).Zero();
	}
}

void Matrix::Reallocate(int _N)
{
	N = _N;
	val.Reallocate(N);
	for(int i(0); i < N; ++i) {
		val(i).Reallocate(N);
	}
}

void Matrix::Resize(int _N)
{
	N = _N;
	val.Resize(N);
	for(int i(0); i < N; ++i) {
		val(i).Resize(N);
	}
}

double& Matrix::operator()(int i, int j)
{
	return val(i)(j);
}

const double& Matrix::operator()(int i, int j) const
{
	return val(i)(j);
}

Matrix& Matrix::operator=(const Matrix& m)
{
	N = m.GetM();
	for(int i(0); i < N; ++i) {
		val(i) = m.val(i);
	}
	return *this;
}

ostream& operator<<(ostream& out, const Matrix& m)
{
	int N(m.GetM());
	for(int i(0), j; i < N; ++i) {
		for(j = 0; j < N; ++j) {
			out << m(i,j) << "\t";
		}
		out << endl;
	}
	return out;
}

#endif
