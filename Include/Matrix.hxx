
#ifndef MATRIX_HXX
#define MATRIX_HXX

#include "define.hxx"
#include <vector>
#include <iostream>

using namespace std;


class Matrix
{
	public:
		Matrix();
		Matrix(const Matrix& m);
		Matrix(int _N, int _M);
		Matrix(int _N);
		int GetN() const; // N number of ligne
		int GetM() const; // M number of row
		void Zero();
		void Reallocate(int _N);
		void Reallocate(int _N,int _M);
		void Clear();
		precision& operator()(int i, int j);
		const precision& operator()(int i, int j) const;
		Matrix& operator=(const Matrix& m);
		Matrix operator+(const Matrix&);
	protected:
		int N,M;
		vector<vector<precision> > val;
};

ostream& operator<<(ostream& out, const Matrix& m);
Matrix operator*(precision,const Matrix&);



#endif