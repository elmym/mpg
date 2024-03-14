#include "first.h"

double Solver::filtration_velocity() {
	return _k * _delta_p / _mu / _L;
}

double Solver::true_velocity() {
	double u = filtration_velocity();
	return u /_m;
}

double Solver:: time_between_galeries() {
	return _m * _mu * _L * _L / _k / (_delta_p);
}

double Reservoir_Solver::time_between(double m, double L, double mu, double k, double delta_p) {
	double sum = 0;
	for (int i = 0; i < _grid_size - 1; i++) {
		sum += (_grid[i + 1] - _grid[i]) * (_grid[i + 1] - _grid[i]) / (-_p[i + 1] + _p[i]);
	}
	double u = ((delta_p) / sum) / mu / L / m;
	return L / u / k;
}

static void thomas_algorithm(const std::vector<double>& a,
	const std::vector<double>& b,
	const std::vector<double>& c,
	const std::vector<double>& d,
	std::vector<double>& f) {
	int n = b.size();
	std::vector<double> v(n);
	std::vector<double> u(n);
	// forward stroke
	v[0] = -c[0] / b[0];
	u[0] = d[0] / b[0];
	for (int i = 1; i < n - 1; i++) {
		v[i] = -c[i] / (b[i] + a[i] * v[i - 1]);
		u[i] = (a[i] * u[i - 1] - d[i]) / (-b[i] - a[i] * v[i - 1]);
	}
	v[n - 1] = 0;
	u[n - 1] = -(a[n - 1] * u[n - 2] - d[n - 1]) / (b[n - 1] + a[n - 1] * v[n - 2]);
	//reversal
	f[n - 1] = u[n - 1];
	for (int i = n - 1; i > 0; i--) {
		f[i - 1] = v[i - 1] * f[i] + u[i - 1];
	}
}

Reservoir_Solver::Reservoir_Solver(std::string input, double a, double k1, double k2) :_a(a), _k1(k1), _k2(k2) {
	std::ifstream fs(input);
	copy(std::istream_iterator<double>(fs), std::istream_iterator<double>(), back_inserter(_grid));
	_grid_size = _grid.size();
	fs.close();
}

Reservoir_Solver::Reservoir_Solver(size_t N, double L, double a, double k1, double k2) :_a(a), _k1(k1), _k2(k2) {
	_grid_size = N + 1;
	double len = 0;
	double h = L / N;
	for (int i = 0; i < N; i++) {
		_grid.push_back(len);
		len += h;
	}
	_grid.push_back(len);
}
double Reservoir_Solver::exact_p(double x) {
	if (_a == 1) {
		double p = 1 - x;
		return p;
	}
	else {
		double pa = _k1 / (_k1 + _k2);
		if (x < _a) 
			return 1 - (1 - pa) * x / _a;
		else if (x == _a) 
			return pa;
		else 
			return pa - (pa) * (x - _a) / (1 - _a);
	}
}

void Reservoir_Solver::find_p() {
	_p.resize(_grid_size);
	std::vector<double> A(_grid_size);
	std::vector<double> B(_grid_size);
	std::vector<double> C(_grid_size);
	A[0] = 0; B[0] = 1; C[0] = 0;
	if (_a == 1) {
		for (int i = 1; i < _grid_size - 1; i++) {
			A[i] = -_grid[i] + _grid[i + 1];
			B[i] = -_grid[i + 1] + _grid[i - 1];
			C[i] = -_grid[i - 1] + _grid[i];
		}
	}
	else {
		int j = 1;
		for (int i = 0; i < _grid_size; i++) {
			if (_grid[i] > _a) {
				j = i;
				break;
			}
		}
		for (int i = 1; i < j; i++) {
			A[i] = _k1 * (_grid[i] - _grid[i - 1]) / _grid[i] * _a;
			C[i] = _k1 * (_grid[i + 1] - _grid[i]) / _grid[i] * _a;
			B[i] = -A[i] - C[i];
		}
		for (int i = j; i < _grid_size - 1; i++) {
			A[i] = _k2 * (_grid[i] - _grid[i - 1]) / (_grid[i] - _a) * (1 - _a);
			C[i] = _k2 * (_grid[i + 1] - _grid[i]) / (_grid[i] - _a) * (1 - _a);
			B[i] = -A[i] - C[i];
		}
	}
	A[_grid_size - 1] = 0;
	B[_grid_size - 1] = 1;
	C[_grid_size - 1] = 0;

	std::vector<double> D(_grid_size);
	D[0] = 1;

	thomas_algorithm(A, B, C, D, _p);
}


double Reservoir_Solver::norm() {
	find_p();
	double norm = 0;
	for (int i = 0; i < _grid_size; i++) {
		double delta = _p[i] - exact_p(_grid[i]);
		norm += delta * delta;
	}
	return std::sqrt(norm / _grid_size);
}
