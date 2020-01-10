//Pang 8.3
//Molecular Dynamics

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <array>
#include <vector>

//Constants
const double kb = 1.38065e-23;	//SI
const double e0 = 8.854e-12;	//SI
const double PI = 3.14159265358979323846;
const double r0 = 0.321e-10;	//meters
const double r_eq = 2.75e-10;	//meters
const double V0 = 1.746e-16;	//Joules

//System parameters
const double T = 300;
const double totTime = 1e-12;
const int stepNum = 100;
const double t = totTime / stepNum;

//Structs defining ions containing mass and lists of r, v, and a as 3D vectors
struct Na {
	const double mass = 3.8175e-26;
	std::vector<std::array<double, 3>> r;
	std::vector<std::array<double, 3>> v;
	std::vector<std::array<double, 3>> a;
};

struct Cl {
	const double mass = 5.8119e-26;
	std::vector<std::array<double, 3>> r;
	std::vector<std::array<double, 3>> v;
	std::vector<std::array<double, 3>> a;
};

//Returns a 3D array containing random vx, vy, and vz according to maxwell distribution
std::array<double,3> maxwell(double m, double T);

//Returns accelleration as a function of position for sodium and chlorine ions
std::array<double, 3> gn(std::vector<Na> const& s, std::vector<Cl> const& c, int p_index);
std::array<double, 3> gc(std::vector<Na> const& s, std::vector<Cl> const& c, int p_index);
int main()
{
	//Create (n^3)/2 sodium ions and (n^3)/2 chlorine ions with initial positions in equilibrium (a[0] = (0,0,0)),
	//and initial velocities according to maxwell distribution
	int n = 4;
	std::vector<Na> Sodium;
	std::vector<Cl> Chlorine;
	std::cout << "Creating particles...\n";
	for (int k = 0; k < n; k++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < n; i++)
			{
				//Create and Assign positions for Sodium and Chlorine ions in 2 interlaced FCC lattices
				if (k % 2 == 0) 
				{
					if (j % 2 == 0)
					{
						if (i % 2 == 0)
						{
							Sodium.push_back(Na());
							Sodium.back().r.push_back({ i * r_eq,j * r_eq,k * r_eq });
						}
						else
						{
							Chlorine.push_back(Cl());
							Chlorine.back().r.push_back({ i * r_eq,j * r_eq,k * r_eq });
						}
					}
					else
					{
						if (i % 2 == 0)
						{
							Chlorine.push_back(Cl());
							Chlorine.back().r.push_back({ i * r_eq,j * r_eq,k * r_eq });
						}
						else
						{
							Sodium.push_back(Na());
							Sodium.back().r.push_back({ i * r_eq,j * r_eq,k * r_eq });
						}
					}
				}
				else 
				{
					if (j % 2 == 0)
					{
						if (i % 2 == 0)
						{
							Chlorine.push_back(Cl());
							Chlorine.back().r.push_back({ i * r_eq,j * r_eq,k * r_eq });
						}
						else
						{
							Sodium.push_back(Na());
							Sodium.back().r.push_back({ i * r_eq,j * r_eq,k * r_eq });
						}
					}
					else
					{
						if (i % 2 == 0)
						{
							Sodium.push_back(Na());
							Sodium.back().r.push_back({ i * r_eq,j * r_eq,k * r_eq });
						}
						else
						{
							Chlorine.push_back(Cl());
							Chlorine.back().r.push_back({ i * r_eq,j * r_eq,k * r_eq });
						}
					}
				}
			}
		}
	}

	int sNum = Sodium.size();
	int cNum = Chlorine.size();


	for (int i = 0; i < sNum; i++)
	{
		//Assign random velocity and acceleration for sodium
		Sodium[i].v.push_back(maxwell(Sodium[i].mass, T));
		Sodium[i].a.push_back({ 0.0,0.0,0.0 });
	}

	for (int i = 0; i < cNum; i++)
	{
		//Assign random velocity and acceleration for chlorine
		Chlorine[i].v.push_back(maxwell(Chlorine[i].mass, T));
		Chlorine[i].a.push_back({ 0.0,0.0,0.0 });
	}

	std::cout << "Implementing Verlet algorithm...\n";
	//Implement Verlet algorithm
	for (int m = 1; m < stepNum; m++)
	{
		//Update r for all particles
		for (int i = 0; i < sNum; i++)
		{
			//Update r[m] with previous r, v, and a.
			double x = Sodium[i].r[m - 1][0] + t * Sodium[i].v[m - 1][0] + t * t / 2.0 * Sodium[i].a[m - 1][0];
			double y = Sodium[i].r[m - 1][1] + t * Sodium[i].v[m - 1][1] + t * t / 2.0 * Sodium[i].a[m - 1][1];
			double z = Sodium[i].r[m - 1][2] + t * Sodium[i].v[m - 1][2] + t * t / 2.0 * Sodium[i].a[m - 1][2];
			Sodium[i].r.push_back({ x,y,z });
		}
		for (int i = 0; i < cNum; i++)
		{
			//Update r[m] with previous r, v, and a.
			double x = Chlorine[i].r[m - 1][0] + t * Chlorine[i].v[m - 1][0] + t * t / 2.0 * Chlorine[i].a[m - 1][0];
			double y = Chlorine[i].r[m - 1][1] + t * Chlorine[i].v[m - 1][1] + t * t / 2.0 * Chlorine[i].a[m - 1][1];
			double z = Chlorine[i].r[m - 1][2] + t * Chlorine[i].v[m - 1][2] + t * t / 2.0 * Chlorine[i].a[m - 1][2];
			Chlorine[i].r.push_back({ x,y,z });
		}
		//Update a with new r
		for (int i = 0; i < sNum; i++)
			Sodium[i].a.push_back(gn(Sodium, Chlorine, i));
		for (int i = 0; i < cNum; i++)
			Chlorine[i].a.push_back(gc(Sodium, Chlorine, i));

		//Update v
		for (int i = 0; i < sNum; i++)
		{
			double x = Sodium[i].v[m - 1][0] + t / 2.0 * (Sodium[i].a[m - 1][0] + Sodium[i].a[m][0]);
			double y = Sodium[i].v[m - 1][1] + t / 2.0 * (Sodium[i].a[m - 1][1] + Sodium[i].a[m][1]);
			double z = Sodium[i].v[m - 1][2] + t / 2.0 * (Sodium[i].a[m - 1][2] + Sodium[i].a[m][2]);
			Sodium[i].v.push_back({ x,y,z });
		}
		for (int i = 0; i < cNum; i++)
		{
			double x = Chlorine[i].v[m - 1][0] + t / 2.0 * (Chlorine[i].a[m - 1][0] + Chlorine[i].a[m][0]);
			double y = Chlorine[i].v[m - 1][1] + t / 2.0 * (Chlorine[i].a[m - 1][1] + Chlorine[i].a[m][1]);
			double z = Chlorine[i].v[m - 1][2] + t / 2.0 * (Chlorine[i].a[m - 1][2] + Chlorine[i].a[m][2]);
			Chlorine[i].v.push_back({ x,y,z });
		}
	}

	std::cout << "Calculating kinetic energies...\n";
	//Calculate Kinetic energy of the crystal from E = (mv^2)/2
	std::vector<double> Ek;
	for (int m = 0; m < stepNum; m++)
	{
		Ek.push_back(0.0);
		for (int i = 0; i < sNum; i++)
		{
			double v = sqrt(pow(Sodium[i].v[m][0], 2) + pow(Sodium[i].v[m][1], 2) + pow(Sodium[i].v[m][2], 2));
			Ek[m] += Sodium[i].mass * v * v / 2;
		}
		for (int i = 0; i < cNum; i++)
		{
			double v = sqrt(pow(Chlorine[i].v[m][0], 2) + pow(Chlorine[i].v[m][1], 2) + pow(Chlorine[i].v[m][2], 2));
			Ek[m] += Chlorine[i].mass * v * v / 2;
		}
	}

	//Output Data

	std::cout << "Number of Na ions: " << sNum << std::endl;
	std::cout << "Number of Cl ions: " << cNum << std::endl;
	std::cout << "Initial kinetic energy: " << Ek.front() << std::endl;
	std::ofstream fout ("Ek.txt");
	fout << "\"[";
	for (int m = 0; m < stepNum; m++)
	{
		fout << Ek[m] << ",";
	}
	//Delete last ',' and add ']'
	long pos = fout.tellp();
	fout.seekp(pos - 1.0);
	fout << "]\"";

	fout.close();
	fout.open("data.txt");
	fout << "\"[";
	for (int m = 0; m < stepNum; m++)
	{
		fout << "[[";
		for (int i = 0; i < sNum - 1; i++)
		{
			fout << "[" << Sodium[i].r[m][0] << "," << Sodium[i].r[m][1] << "," << Sodium[i].r[m][2] << "],";
		}
		fout << "[" << Sodium[sNum - 1].r[m][0] << "," << Sodium[sNum - 1].r[m][1] << "," << Sodium[sNum - 1].r[m][2] << "]],[";
		for (int i = 0; i < cNum - 1; i++)
		{
			fout << "[" << Chlorine[i].r[m][0] << "," << Chlorine[i].r[m][1] << "," << Chlorine[i].r[m][2] << "],";
		}
		fout << "[" << Chlorine[cNum - 1].r[m][0] << "," << Chlorine[cNum - 1].r[m][1] << "," << Chlorine[cNum - 1].r[m][2] << "]]],";
	}
	//Delete last ',' and add ']\"'
	pos = fout.tellp();
	fout.seekp(pos - 1.0);
	fout << "]\"";

	Sodium.clear();
	Chlorine.clear();

	return 0;
}

std::array<double,3> maxwell(double m, double T)
{
	
	//Initialize gaussian generator
	double sigma = sqrt(kb * T / m);
	static std::default_random_engine generator;
	static std::normal_distribution<double> gauss(0, 1);

	std::array<double, 3> v;
	for (int i = 0; i < 3; i++)
		v[i] = sigma * gauss(generator);

	return v;
}

std::array<double, 3> gn(std::vector<Na> const &s, std::vector<Cl> const &c, int p_index)
{
	//Calculate potential for Na  at p_index
	std::array<double,3> V = { 0 };
	int n = s.size();
	int m = c.size();
	double constant = pow(1.602e-19, 2) / 4 / PI / e0;
	std::array<double, 3> pos = s[p_index].r.back();
	//Sum potential from all other Na ions
	for (int i = 0; i < n; i++)
	{
		if (i != p_index)
		{
			double r = sqrt(pow(pos[0] - s[i].r.back()[0], 2) + pow(pos[1] - s[i].r.back()[1], 2) + pow(pos[2] - s[i].r.back()[2], 2));
			for (int d = 0; d < 3; d++)
				V[d] += constant * (pos[d] - s[i].r.back()[d]) / pow(r, 3);
		}
	}

	//Sum potential from all Cl ions
	for (int i = 0; i < m; i++)
	{
		double r = sqrt(pow(pos[0] - c[i].r.back()[0], 2) + pow(pos[1] - c[i].r.back()[1], 2) + pow(pos[2] - c[i].r.back()[2], 2));
		for (int d = 0; d < 3; d++)
			V[d] += -1 * constant * (pos[d] - c[i].r.back()[d]) / pow(r, 3) + V0 * exp(-r / r0) * (pos[d] - c[i].r.back()[d]) / r / r0;
	}

	//Divide by mass to get acceleration
	for (int d = 0; d < 3; d++)
		V[d] /= s[0].mass;

	//Return V
	return V;
}

std::array<double, 3> gc(std::vector<Na> const& s, std::vector<Cl> const& c, int p_index)
{
	//Calculate potential for Cl  at p_index
	std::array<double, 3> V = { 0 };
	int n = s.size();
	int m = c.size();
	double constant = pow(1.602e-19, 2) / 4 / PI / e0;
	std::array<double, 3> pos = c[p_index].r.back();
	//Sum potential from all other Cl ions
	for (int i = 0; i < m; i++)
	{
		if (i != p_index)
		{
			double r = sqrt(pow(pos[0] - c[i].r.back()[0], 2) + pow(pos[1] - c[i].r.back()[1], 2) + pow(pos[2] - c[i].r.back()[2], 2));
			for (int d = 0; d < 3; d++)
				V[d] += constant * (pos[d] - c[i].r.back()[d]) / pow(r, 3);
		}
	}

	//Sum potential from all Na ions
	for (int i = 0; i < m; i++)
	{
		double r = sqrt(pow(pos[0] - s[i].r.back()[0], 2) + pow(pos[1] - s[i].r.back()[1], 2) + pow(pos[2] - s[i].r.back()[2], 2));
		for (int d = 0; d < 3; d++)
			V[d] += -1 * constant * (pos[d] - s[i].r.back()[d]) / pow(r, 3) + V0 * exp(-r / r0) * (pos[d] - s[i].r.back()[d]) / r / r0;
	}
	//Divide by mass to get acceleration
	for (int d = 0; d < 3; d++)
		V[d] /= c[0].mass;

	//Return V
	return V;
}