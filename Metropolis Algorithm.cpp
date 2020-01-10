//Pang 10.5
//Metropolis Algorithm ferromagnetic Ising model

#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <fstream>

//global variables and constants
int M = 0;		//number of steps completed
int accept = 0;		//number of accepted state changes
const double kb = 1.3807e-23; // J/K
const double T = 1; //Kelvin
const double J = 1.21e-21;	//Exchange interaction strength of Iron in Joules
const double B = 0;	//External magnetic field strength
const int N = 100;		//Lattice length and width
const int accept_thresh = 20000;	//threshold of accepted states before program terminates
const int h = 100;		//Number of steps to skip before measuring magnitization

//Latice point class definition
class lat_pt {
private:
	int s;		//spin
	int sj_sum; //Sum of surrounding spins
public:
	//Constructor
	lat_pt();

	//sets sj_set by adding up the s values from the four nearest neighbors
	void sj_set(lat_pt& p1, lat_pt& p2, lat_pt& p3, lat_pt& p4);
	
	//Update the sj_sum of lattice points surrounding a recently flipped point
	void sj_update(int update_sign, lat_pt& p1, lat_pt& p2, lat_pt& p3, lat_pt& p4);

	//Change state of s
	void s_flip();

	//Returns value of s
	int s_get() const;
	
	//Return the change in Hamiltonian for a given lattice point switching spin states
	//Set B to zero if you would like to analyze the behavior in the absence of an external magnetic field
	double delta_H() const;
};

//Random generator used by random_num functions
std::mt19937& generator();

//Return random double in [min,max]
double random_num(double min, double max);

//Return random int in [min,max]
int random_num(int min, int max);

//Return either -1 or 1 with 50% probability
int random_spin();

//Implements the metropolis algorithm by updating the spin configuration for a single site
//n1, n2, n3, and n4 are the four nearest neigbors that must be passed in order
//to update their sj_sum if the spin of lattice point pt flips.
void metropolis(lat_pt& pt, lat_pt& n1, lat_pt& n2, lat_pt& n3, lat_pt& n4);

int main() {


	//Build lattice with 2D std::vector
	std::vector<std::vector<lat_pt>> lattice;
	for (int x = 0; x < N; x++)
	{
		lattice.push_back(std::vector<lat_pt>());
		for (int y = 0; y < N; y++)
		{
			lattice[x].push_back(lat_pt());
		}
	}

	//set sj_sum for all lattice points besides boundaries
	for (int x = 1; x < N-1; x++)
	{
		for (int y = 1; y < N-1; y++)
		{
			lattice[x][y].sj_set(lattice[x + 1][y], lattice[x - 1][y], lattice[x][y + 1], lattice[x][y - 1]);
		}
	}

	//set boundaries according to periodic boundary conditions
	//set sj_sum for boundaries besides corners
	for (int i = 1; i < N - i; i++)
	{
		lattice[i][0].sj_set(lattice[i + 1][0], lattice[i - 1][0], lattice[i][1], lattice[i][N - 1]);
		lattice[i][N - 1].sj_set(lattice[i + 1][N - 1], lattice[i - 1][N - 1], lattice[i][0], lattice[i][N - 2]);
		lattice[0][i].sj_set(lattice[0][i + 1], lattice[0][i - 1], lattice[1][i], lattice[N - 1][i]);
		lattice[N - 1][i].sj_set(lattice[N - 1][i + 1], lattice[N - 1][i - 1], lattice[0][i], lattice[N - 2][i]);

	}

	//set sj_sum for corners
	lattice[0][0].sj_set(lattice[1][0], lattice[0][1], lattice[0][N - 1], lattice[N - 1][0]);
	lattice[N - 1][0].sj_set(lattice[N - 1][1], lattice[N - 2][0], lattice[N - 1][N - 1], lattice[0][0]);
	lattice[0][N - 1].sj_set(lattice[1][N - 1], lattice[0][N - 2], lattice[N - 1][N - 1], lattice[0][0]);
	lattice[N - 1][N - 1].sj_set(lattice[N - 2][N - 1], lattice[N - 1][N - 2], lattice[N - 1][0], lattice[0][N - 1]);



	//Implement Metropolis algorithm, changing lattice points randomly
	double s_sigma = 0; //sum of all s values
	int conf_num = 0; //number of configurations used to calculate magnetization
	double magnetization;

	//run metropolis algorithm until accept > accept_thresh
	while (accept <= accept_thresh)
	{
		for (int i = 0; i < h; i++)
		{
			int x = random_num(0, N - 1);
			int y = random_num(0, N - 1);
			if (x != 0 && x != N - 1 && y != 0 && y != N - 1)
				metropolis(lattice[x][y], lattice[x + 1][y], lattice[x - 1][y], lattice[x][y + 1], lattice[x][y - 1]);
			else if (x == 0 && y != 0 && y != N - 1)
				metropolis(lattice[x][y], lattice[x + 1][y], lattice[N - 1][y], lattice[x][y + 1], lattice[x][y - 1]);
			else if (x == N - 1 && y != 0 && y != N - 1)
				metropolis(lattice[x][y], lattice[0][y], lattice[N - 1][y], lattice[x][y + 1], lattice[x][y - 1]);
			else if (x != 0 && x != N - 1 && y == 0)
				metropolis(lattice[x][y], lattice[x + 1][y], lattice[x - 1][y], lattice[x][y + 1], lattice[x][N - 1]);
			else if (x != 0 && x != N - 1 && y == N - 1)
				metropolis(lattice[x][y], lattice[x + 1][y], lattice[x - 1][y], lattice[x][0], lattice[x][y - 1]);
			else if (x == 0 && y == 0)
				metropolis(lattice[0][0], lattice[1][0], lattice[0][1], lattice[0][N - 1], lattice[N - 1][0]);
			else if (x == N - 1 && y == 0)
				metropolis(lattice[N - 1][0], lattice[N - 1][1], lattice[N - 2][0], lattice[N - 1][N - 1], lattice[0][0]);
			else if (x == 0 && y == N - 1)
				metropolis(lattice[0][N - 1], lattice[1][N - 1], lattice[0][N - 2], lattice[N - 1][N - 1], lattice[0][0]);
			else
				lattice[N - 1][N - 1].sj_set(lattice[N - 2][N - 1], lattice[N - 1][N - 2], lattice[N - 1][0], lattice[0][N - 1]);


		}
	
		//After h lattice points have attempted to flip, sum all s values
		for (int x = 0; x < N; x++)
		{
			for (int y = 0; y < N; y++)
			{
				s_sigma += lattice[x][y].s_get();
			}
		}
		conf_num++;
	}

	magnetization = s_sigma / (double(N) * N) / conf_num;


	double accept_percent = 100.0 * accept / M;
	std::cout << "Temperature: " << T << " Kelvin\n";
	std::cout << "Exchange interaction strength: " << J << std::endl;
	std::cout << "Magnetization: " << magnetization << std::endl;
	std::cout << "Number of configurations measured: " << conf_num << std::endl;
	std::cout << "Total number of accepted state changes: " << accept << std::endl;
	std::cout << "Total number of steps: " << M << std::endl;
	std::cout << "Acceptance rate = " << accept_percent << "%\n\n";

	//Output final configuration in JSON format to be read by Mathematica ListPlot
	std::ofstream fout("data.txt");
	fout << "\"[[";
	for (int x = 0; x < N; x++)
	{
		for (int y = 0; y < N; y++)
		{
			if (lattice[x][y].s_get() > 0)
				fout << "[" << x << "," << y << "],";
		}
	}
	long pos = fout.tellp();
	fout.seekp(pos - 1.0);
	fout << "],[";
	for (int x = 0; x < N; x++)
	{
		for (int y = 0; y < N; y++)
		{
			if (lattice[x][y].s_get() < 0)
				fout << "[" << x << "," << y << "],";
		}
	}
	pos = fout.tellp();
	fout.seekp(pos - 1.0);
	fout << "]]\"";
	fout.close();
	return 0;
}





int random_spin() {
	int spin = random_num(0, 1);
	if (spin)
		return spin;
	else
		return -1;
}

 std::mt19937& generator() {
	static std::random_device rd;
	static std::mt19937 gen(rd());
	return gen;
}

double random_num(double min, double max) {
	static std::uniform_real_distribution<double> dist{};
	using parm_t = decltype(dist)::param_type;
	return dist(generator(), parm_t{ min, std::nextafter(max, DBL_MAX) });
}

int random_num(int min, int max) {
	static std::uniform_int_distribution<int> dist{};
	using parm_t = decltype(dist)::param_type;
	return dist(generator(), parm_t{ min, max });
}


void metropolis(lat_pt& pt, lat_pt& n1,lat_pt& n2,lat_pt& n3,lat_pt& n4) {
	
	//flip s and calculate probability of change of state
	pt.s_flip();
	double q = exp(-pt.delta_H() / kb / T);

	if (q >= random_num(0.0, 1.0))
	{
		//accept change of state
		//update the sj_sum value of 4 nearest neighbors
		pt.sj_update(pt.s_get(),n1, n2, n3, n4);
		accept++;
	}
	else
	{
		//reject change of state
		pt.s_flip();
	}

	//increment M;
	M++;
}


//Methods for lat_pt class
lat_pt::lat_pt() {
	s = random_spin();
}

void lat_pt::sj_set(lat_pt& p1, lat_pt& p2, lat_pt& p3, lat_pt& p4) {
	sj_sum = p1.s + p2.s + p3.s + p4.s;
}

void lat_pt::sj_update(int update_sign, lat_pt& p1, lat_pt& p2, lat_pt& p3, lat_pt& p4) {
	p1.sj_sum += 2 * update_sign;
	p2.sj_sum += 2 * update_sign;
	p3.sj_sum += 2 * update_sign;
	p4.sj_sum += 2 * update_sign;
}

void lat_pt::s_flip()
{
	s *= -1;
}

int lat_pt::s_get() const
{
	return s;
}

double lat_pt::delta_H() const
{
	return -2 * J * s * sj_sum - 2 * B * s;
}
