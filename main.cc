#include "linalg/Linalg.hxx"


#include "velocity.cxx"

using namespace linalg;


int main()
{

	Velocity v(10, 10, 5.0, 5.0);
	v.ChampsCirculaire(2.5, 2.5, 0.5);
	v.WriteGnuPlot("velocity.dat");
	// plot "velocity.dat" u 1:2:3:4 w vec
	cout << "ok cela tourne!";


	return 0;
}
