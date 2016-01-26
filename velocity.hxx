#ifndef FILE_VELOCITY_HXX


namespace linalg
{
class Velocity

{
private:
	Vector<Vector<double> > VX,VY;
	int Nx,Ny;
	double L,H,Delta_x,Delta_y;

public:
	Velocity(int Nx,int Ny,double L,double H);
	void ChampsCirculaire(double Xcenter,double Ycenter, double intensite);
	void WriteGnuPlot(const string& nom);
	double& GetVX(int i, int j);
	double& GetVY(int i, int j);

};
}
#define FILE_VELOCITY_HXX
#endif
