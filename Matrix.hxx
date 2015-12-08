
#ifndef MATRIX_HXX
#define MATRIX_HXX

#include "linalg/Linalg.hxx"
#include <iostream>

using namespace std;
using namespace linalg;

//coucou julien

class Matrix
{
	public:
		Matrix();
		Matrix(const Matrix& m);
		Matrix(int _N);
		int GetM() const;
		void Zero();
		void Reallocate(int _N);
		void Resize(int _N);
		double& operator()(int i, int j);
		const double& operator()(int i, int j) const;
		Matrix& operator=(const Matrix& m);
	protected:
		int N;
		Vector<Vector<double> > val;
};

ostream& operator<<(ostream& out, const Matrix& m);


#endif
