//Pang 11.2
//Binary Genetic algorithm

#include <iostream>
#include <random>
#include <cmath>
#include <vector>


const int prec = 16;	//precision - number of bits
const int param = 2;	//number of parameters in optimized function
const double PI = 3.14159265358979323846;
const double mutatePercent = 0.05;	//percent of bits that are mutated each generation
const int mutateSkip = 10;	//number of chromosomes to avoid mutating
const int generations = 150;	//number of generations to produce
const int population = 500;		//starting population

//Random generator used by random_num functions
std::mt19937& generator();

//Return random double in [min,max]
double random_num(double min, double max);

//Return random int in [min,max]
int random_num(int min, int max);

//cost function
double f(double x, double y) { return y * sin(4 * PI * x) + 2 * x * cos(8 * PI * y); }

//class representing a chromosome
class chrom {
private:
	//array containing encoded parameters
	bool w[param][prec];

	//value of function at given parameter space
	double cost;

public:
	//constructor randomly assigns bits to w
	chrom();

	//encodes an array of parameters as binary in w
	void encode(const double p[param]);

	//decodes w parameters to array of doubles
	double* decode() const;

	//perform single point crossover and produce 2 offspring
	std::vector<chrom> reproduce(chrom& mate);

	//update w for reproduce method
	void change_w(int i, int j, bool bit);

	//mutate bits randomly
	void mutate();

	//sets cost
	void set_cost();

	//returns cost
	double get_cost() const;

	//overide < operator for sort function
	bool operator < (const chrom& a) const
	{
		return (cost < a.get_cost());
	}

	void display_chrom() const;
};

int main()
{
	std::vector<chrom> parents(population);		//vector to hold parents
	std::vector<chrom> children(population);	//vector to hold children

	//implement binary genetic algorithm
	for (int g = 0; g < generations; g++)
	{
		//shuffle gene pool
		std::shuffle(parents.begin(), parents.end(), generator());
		//compare consecutive chromosomes and delete the one with highest cost
		for (int i = 0; i < population / 2; i++)
		{
			if (parents[i] < parents[i + 1])
			{
				parents.erase(parents.begin() + i + 1);
			}
			else parents.erase(parents.begin() + i);
		}

		

		//keep the 4 best parents in the next generation and reproduce
		sort(parents.begin(), parents.end());
		for (int i = 0; i < 4; i++)
			children[i] = parents[i];
		std::shuffle(parents.begin(), parents.end(), generator());
		for (int i = 4; i < population / 2; i += 2)
		{
			swap_ranges(children.begin() + i, children.begin() + i + 1, parents[i].reproduce(parents[i + 1]).begin());
		}
		//shuffle and repeat to restore original population
		std::shuffle(parents.begin(), parents.end(), generator());
		//add new children to end of vector to prevent excessive memory reallocation
		for (int i = 0; i < population / 2; i += 2)
		{
			swap_ranges(children.end() - i - 1, children.end() - i, parents[i].reproduce(parents[i + 1]).begin());
		}


		//sort children to prepare for mutation
		sort(children.begin(), children.end());
		//mutate top configurations, but undo if cost increases
		for (int i = 0; i < mutateSkip; i++)
		{
			chrom temp = children[i];
			children[i].mutate();
			if (temp < children[i])
				children[i] = temp;
		}
		//mutate rest of children
		for (int i = mutateSkip; i < population; i++)
			children[i].mutate();


		//prepare for next generation
		parents = children;
		if (!((g + 1) % 25))
		{
			std::cout << "Completed generation " << g + 1 << std::endl;
			parents[0].display_chrom();
		}
	}

	//display top 5% chromosomes
	std::cout << "Top 5% chromosomes:\n-------------------------------\n";
	sort(parents.begin(), parents.end());
	for (int i = 0; i < population*.05; i++)
		parents[i].display_chrom();
	std::cout << "-------------------------------\n\n";

	//Output results
	double* p = parents[0].decode();
	double x = 2 * p[0] - 1;
	double y = 2 * p[1] - 1;
	delete[] p;
	double min = parents[0].get_cost();

	using namespace std;
	cout << "Minumum of " << min << " was found at (" << x << ", " << y << ")" << endl;
	cout << "after " << generations << " generations of " << population << " chromosomes.\n\n";

	return 0;

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



chrom::chrom()
{
	for (int i = 0; i < param; i++)
	{
		for (int j = 0; j < prec; j++)
		{
			w[i][j] = bool(random_num(0, 1));
		}
	}
	set_cost();
}


void chrom::encode(const double p[param])
{
	for (int i = 0; i < param; i++)
	{
		w[i][0] = round(p[i]);
		for (int j = 1; j < prec; j++)
		{
			double sum = 0.0;
			for (int k = 0; k < j; k++)
			{
				sum += pow(2, j - k - 1) * double(w[i][k]);
			}
			w[i][j] = round(pow(2, j)*p[i] - sum);
		}
	}
	//set cost of new configuration
	set_cost();
}

double* chrom::decode() const
{
	double* arr = new double[param];
	
	for (int i = 0; i < param; i++)
	{
		arr[i] = 0.0;
		for (int j = 0; j < prec; j++)
			if (w[i][j])
				arr[i] += 1.0 / pow(2, j + 1);
	}
	return arr;
}

std::vector<chrom> chrom::reproduce(chrom& mate)
{
	std::vector<chrom> c(2);	//2 children
	//let crossover occur somewhere randomly in the chromosome
	int x = random_num(1, prec * param);
	for (int i = 0; i < x; i++)
	{
		c[0].change_w(i / prec, i % prec, w[i / prec][i % prec]);
		c[1].change_w(i / prec, i % prec, mate.w[i / prec][i % prec]);
	}
	for (int i = x; i < prec * param; i++)
	{
		c[1].change_w(i / prec, i % prec, w[i / prec][i % prec]);
		c[0].change_w(i / prec, i % prec, mate.w[i / prec][i % prec]);
	}

	c[0].set_cost();
	c[1].set_cost();
	return c;
}


void chrom::set_cost()
{
	double* p = decode();
	//scale x and y for domain [-1,1]
	double x = 2 * p[0] - 1;
	double y = 2 * p[1] - 1;

	cost = f(x, y);

	delete[] p;
}

double chrom::get_cost() const
{
	return cost;
}

void chrom::change_w(int i, int j, bool bit)
{
	w[i][j] = bit;
}

void chrom::mutate() 
{
	for (int i = 0; i < param; i++)
	{
		for (int j = 0; j < prec; j++)
		{
			if (random_num(0.0, 1.0) < mutatePercent)
				w[i][j] = !w[i][j];
		}
	}
	set_cost();
}

void chrom::display_chrom() const
{
	double* p = decode();
	//scale x and y for domain [-1,1]
	double x = 2 * p[0] - 1;
	double y = 2 * p[1] - 1;
	std::cout << "Location: (" << x << ", " << y << ")\n";
	delete[] p;
	std::cout << "Cost: " << cost << "\n\n";
}
