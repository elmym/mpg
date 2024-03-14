#include "first.h"



int main()
{
	double L = 100;
	double k = 1e-12;
	double m = 0.2;
	double delta_p = 1e6;
	double mu = 1e-3;
	double p_0 = 1;
	double p_1 = 0;

	Reservoir_Solver worker1("input.txt", 1, 1, 1);
	Reservoir_Solver worker2(5, 1, 1, 1, 1);

	std::cout << worker1.norm() << std::endl;
	std::cout << worker2.norm() << std::endl;

	Solver s(k, mu, m, L, delta_p);
	double u = s.filtration_velocity();
	double v = s.true_velocity();
	double t = s.time_between_galeries();
	double tau = worker2.time_between(m, L, mu, k, delta_p);

	std::cout << "filtration velocity= " << u << ", true velocity= " << v << std::endl;
	std::cout <<"time between galeries exact= " << t/86400 << " days, time betweem galeries numeric= "<<tau/86400<<" days."<<std::endl;

	Reservoir_Solver worker3(10, 1, 0.5, 1.0, 0.1);
	std::cout << worker3.norm() << std::endl;
}
