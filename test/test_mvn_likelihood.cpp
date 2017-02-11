#include <iostream>
#include <fstream>
#include <random>

#include <statistics/multivariatenormal.h>
#include <np_suffies.h>

using namespace std;
using namespace Eigen;

/**
 * Test the multivariate normal distribution. The test uses a predefined mean and covariance matrix and checks with
 * particular data items if the calculated probability is the same as the one calculated beforehand using octave.
 */
int main() {
	default_random_engine generator;

	Vector2d mean(1,1);
	Matrix2d covar;
	covar << 2,0,1,2;

	Suffies_MultivariateNormal suffies(mean.size());
	suffies.mu = mean;
	suffies.sigma = covar;

	multivariate_normal_distribution distribution(suffies);

	data_t data0 = {1,2};

	double prob = distribution.probability(data0);

	cout << "Probability: " << prob << endl;
	assert (prob - 0.061974 < 0.00001);
	
	data_t data1 = {1,2};

	dataset_t dataset;
	dataset.push_back(&data0);
	dataset.push_back(&data1);

	double prob_dataset = distribution.probability(dataset);

	cout << "Probability: " << prob_dataset << endl;
	assert (prob_dataset - 0.0038409 < 0.00001);
}
