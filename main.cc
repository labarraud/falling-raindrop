#include "linalg/Linalg.hxx"


#include "velocity.cxx"

using namespace linalg;


int main()
{
	Velocity v(50, 50, 5.0, 5.0);
	v.ChampsCirculaire(0.0, 0.0, 0.5);
	v.WriteGnuPlot("velocity.dat");
	// plot "velocity.dat" u 1:2:3:4 w vec
	return 0;
}
