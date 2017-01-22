#include <statistics/multivariatenormal_rnd.h>
#include <iostream>
#include <fstream>
#include <random>

using namespace std;
using namespace Eigen;

int main() {
	//static std::mt19937 generator{ std::random_device{}() };
	default_random_engine generator;

	Vector2d mean(1,1);
	Matrix2d covar;
	covar << 1,0,0.2,1;
			
	multivariate_normal_distribution distribution(mean, covar);

	const int nRolls = 10000;
	ofstream ofile;
	ofile.open("test_multivariate_normal_distribution.data");
	for (int i = 0; i < nRolls; ++i) {
		ofile << distribution(generator).transpose() << endl;
	}
	ofile.close();
}
