/*
 * NavierStokes.cxx
 *
 *  Created on: 4 avr. 2016
 *      Author: tlariviere
 */

#include "../Include/NavierStokes.hxx"


NavierStokes::NavierStokes(const Velocity& _v, const Matrix& _p)
	:	v(_v), p(_p)
{

}


void NavierStokes::WriteVtk(const string& nom) const
{

}

void SolveLaplacian(const Matrix& b)
{
// construire A
	// gradient conjugué...
}
