#ifndef FIRST_HPP
#define FIRST_HPP
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

class Reservoir_Solver {
public:
	Reservoir_Solver(size_t N, double L, double a, double k1, double k2);
	Reservoir_Solver(std::string input, double a, double k1, double k2);
	double exact_p(double x);
	void find_p();
	double norm();
	double time_between(double m, double L, double mu, double k, double delta_p);

private:
	std::vector<double> _grid;
	int _grid_size;
	std::vector<double> _p;
	double _a;
	double _k1;
	double _k2;
};

class Solver {
public:
	Solver(double k, double mu, double m, double L, double delta_p) :_k(k), _mu(mu), _m(m), _L(L), _delta_p(delta_p){}
	double filtration_velocity();
	double true_velocity();
	double time_between_galeries();
private:
	double _k;
	double _mu;
	double _m;
	double _L;
	double _delta_p;
};
#endif
