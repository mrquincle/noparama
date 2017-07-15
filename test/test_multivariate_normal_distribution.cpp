#include <iostream>
#include <fstream>
#include <random>

#include <statistics/multivariatenormal.h>
#include <np_suffies.h>

using namespace std;
using namespace Eigen;

int main() {
	default_random_engine generator;
	generator.seed(23849398);

	Vector2d mean(1,1);
	Matrix2d covar;
	covar << 1,0,0.2,1;
			
	Suffies_MultivariateNormal suffies(mean.size());
	suffies.mu = mean;
	suffies.sigma = covar;

	multivariate_normal_distribution distribution(suffies);

	const int nRolls = 10000;
	ofstream ofile;
	ofile.open("test_multivariate_normal_distribution.data");
	for (int i = 0; i < nRolls; ++i) {
		auto res = distribution(generator);
		ofile << res->mu.transpose() << endl;
	}
	ofile.close();
	std::cout << "Wrote data to test_multivariate_normal_distribution.data" << std::endl;
}
