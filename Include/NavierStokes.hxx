/*
 * NavierStokes.hxx
 *
 *  Created on: 4 avr. 2016
 *      Author: tlariviere
 */

#ifndef INCLUDE_NAVIERSTOKES_HXX_
#define INCLUDE_NAVIERSTOKES_HXX_

#include "Matrix.hxx"
#include "velocity.hxx"
#include <iostream>
#include <fstream>
#include <string>

class NavierStokes
{
public:
	NavierStokes(const Velocity& _v, const Matrix& _p);
	void WriteVtk(const string& nom) const;
	void SolveLaplacian(const Matrix& b);
	void Advance(int n, double tn);
private:
	Velocity v;
	Matrix p; // pression
};


#endif /* INCLUDE_NAVIERSTOKES_HXX_ */
