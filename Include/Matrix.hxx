
#ifndef MATRIX_HXX
#define MATRIX_HXX

#include "define.hxx"
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;


enum BoundaryCondition
{
	periodic,
    dirichlet,
    neumann,
};

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
		const BoundaryCondition getbcX0()const;
		const BoundaryCondition getbcXn()const;
		const BoundaryCondition getbcY0()const;
		const BoundaryCondition getbcYn()const;
		const precision getX0()const;
		const precision getXn()const;
		const precision getY0()const;
		const precision getYn() const;
		void SetBoundaryCondition(BoundaryCondition bcX0,precision X0,BoundaryCondition bcXn,precision Xn, BoundaryCondition bcY0, precision Y0,BoundaryCondition bcYn, precision Yn);
		precision norme2();
		precision& operator()(int i, int j);
		precision& operator()(int n);
		const precision operator()(int i, int j) const;
		const precision operator()(int n) const;
		const precision bottom(int i, int j) const;
		const precision top(int i, int j) const;
		const precision right(int i, int j) const;
		const precision left(int i, int j) const;
		Matrix& operator=(const Matrix& m);
		Matrix operator+(const Matrix&);
		precision distnorme2(const Matrix&);
		void Mat2Vec(vector<precision>& out) const;
		void Vec2Mat(const vector<precision>& in);
		void WriteGnuPlot(const string& nom) const;
		void WriteVtk(const string& nom, precision dx, precision dy) const;

	protected:
		int N,M;
		vector<vector<precision> > val;
		precision X0,Xn,Y0,Yn;
		BoundaryCondition bcX0,bcXn,bcY0,bcYn;

};

ostream& operator<<(ostream& out, const Matrix& m);
Matrix operator*(precision,const Matrix&);
precision VecNorme(const vector<precision>& v);
precision VecNorme(const Matrix& v);
precision DotProduct(const vector<precision>& v1, const vector<precision>& v2);
precision DotProduct(const Matrix& v1, const vector<precision>& v2);
unsigned int Bij(int i, int j, int M);



#endif
