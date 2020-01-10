//Pang 4.3
#include <iostream>
#include <cstdlib>
#include <cmath>

//Declare time steps
const int n = 1000;
const int m = 10;
const double totTime = 4 * pow(10, 7);
const double dt = totTime / n;
//Initialize constants
const double Me = 5.972 * pow(10, 24);
const double Ms = 1.989 * pow(10, 30);
const double M = Me + Ms;
const double G = 6.674 * pow(10, -11);
const double PI = 3.14159265358979323846;
const double e = 0.01671;		//Eccentricity of Earth's orbit
const double a = 1.496 * pow(10, 11);	//Semimajor axis of Earth's orbit

//Set initial conditions of Earth's orbit
const double r0 = 1.471 * pow(10, 11);
const double vr0 = 0;
const double Theta0 = 0;
const double Omega0 = 30300 / (1.471 * pow(10, 11));

//Implement the 3/8-rule Runge-Kutta method
//given initial conditions y0[] and length of y0
double* RK4(double y0[], int len, double t, double dt);

//Returns an array of velocities that depend on position variables and time
//The RK4 function can be applied to any system where y' = g(y,t) as long
//as g is defined accordingly.
//g accepts length of y as a parameter to ensure that memory is not corrupted in
//case of improper definition.
double* g(double y[], int len, double t);


int main()
{
	const int len = 4;
	double r[n + 1];
	double vr[n + 1];
	double Theta[n + 1];
	double Omega[n + 1];
	double y[len];

	//Set up initial conditions
	y[0] = r[0] = r0;
	y[1] = vr[0] = vr0;
	y[2] = Theta[0] = Theta0;
	y[3] = Omega[0] = Omega0;

	//Perform Runge-Kutta integration
	for (int i = 0; i < n; i++)
	{
		double t = dt * i;
		double* temp = RK4(y, len, t, dt);
		for (int j = 0; j < len; j++)
			y[j] = temp[j];
		delete[] temp;
		r[i + 1] = y[0];
		vr[i + 1] = y[1];
		Theta[i + 1] = y[2];
		Omega[i + 1] = y[3];

		//Bring Theta back to region [-pi, pi]
		int np = (int)(Theta[i + 1] / (2 * PI) + 0.5);
		Theta[i + 1] -= 2 * PI * np;
	}
	double exact;
	//Compare results to exact values at every m steps;
	using namespace std;
	for (int i = 0; i < n + 1; i += m)
	{
		cout << scientific;
		exact = a * (1 - e * e) / (1 + e * cos(Theta[i]));
		cout << "Calculated(r,theta): (" << r[i] << ", " << fixed << Theta[i] << scientific << ")\n";
		cout << "Exact(r,theta): (" << exact << ", " << fixed << Theta[i] << ")\n";
		cout << "Percent Error: %" << abs(exact - r[i]) / exact * 100.0 << endl << endl;
	}

	return 0;
}



double* RK4(double y0[], int len, double t, double dt)
{
	//Create array to hold k1-k4 for each parameter in y0
	double* k1;
	double* k2;
	double* k3;
	double* k4;
	double* temp = new double[len];
	//Create k1-k4 for each parameter of y0
	k1 = g(y0, len, t);
	for (int i = 0; i < len; i++)
		temp[i] = y0[i] + dt * k1[i] / 3;
	k2 = g(temp, len, t + dt / 3);
	for (int i = 0; i < len; i++)
		temp[i] = y0[i] - dt * k1[i] / 3 + dt * k2[i];
	k3 = g(temp, len, t + dt * 2 / 3);
	for (int i = 0; i < len; i++)
		temp[i] = y0[i] + dt * k1[i] - dt * k2[i] + dt * k3[i];
	k4 = g(temp, len, t + dt);
	//Generate y[t + dt]
	for (int i = 0; i < len; i++)
		temp[i] = y0[i] + dt / 8 * (k1[i] + 3 * (k2[i] + k3[i]) + k4[i]);
	//Clear k1-k4
	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;
	//return temp; temp must be deallocated after function call in main()
	return temp;
}

double* g(double y[], int len, double t)
{
	//y = [r, vr, Theta, Omega]
	// dr/dt = vr
	// dvr/dt = -G(Me + Ms)/r^2 + r*Omega^2
	// dTheta/dt = Omega
	// dOmega/dt = -2*vr*Omega/r
	if (4 == len)
	{
		double* v = new double[4];
		v[0] = y[1];
		v[1] = -G * M / (y[0] * y[0]) + y[0] * y[3] * y[3];
		v[2] = y[3];
		v[3] = -2 * y[1] * y[3] / y[0];
		return v;
	}
	std::cout << "Improper definition of g...\n";
	exit(EXIT_FAILURE);
}