#include "../include/GaussianApproximation.h"

#include <iostream>
int main1(int argc, char* argv[]) {

	GaussianApproximation ga(0.00000000000001);
	double r = ga.GetChannelErrorProbability(0, 8);
	
	double t = ga.phi(10.1);
	std::cout << t << std::endl;
	return 0;
}