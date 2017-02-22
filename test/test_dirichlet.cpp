#include <iostream>
#include <fstream>
#include <random>

#include <experimental/filesystem>

#include <np_suffies.h>
#include <statistics/normalinvwishart.h>
#include <statistics/dirichlet.h>

namespace fs = std::experimental::filesystem;

using namespace std;
using namespace Eigen;

/**
 * Test the multivariate normal distribution. The test uses a predefined mean and covariance matrix and checks with
 * particular data items if the calculated probability is the same as the one calculated beforehand using octave.
 */
int main() {
	int _verbosity = 1;
	
	// number of samples to generate
	int N = 100;

	default_random_engine generator;
	
	fout << "Dirichlet distribution" << endl;
	Suffies_Dirichlet suffies_dirichlet;
	suffies_dirichlet.alpha = 1;

	fout << "Normal Inverse Wishart distribution" << endl;
	Suffies_NormalInvWishart suffies_niw(1);
	suffies_niw.mu << 6;
	suffies_niw.kappa = 1.0;
	suffies_niw.nu = 4;
	suffies_niw.Lambda << 0.01;

	normal_inverse_wishart_distribution prior(suffies_niw);

	dirichlet_distribution hyper(suffies_dirichlet, prior);

	bool success;

	const string & path = "output/test";
	success = fs::create_directories(path);
	if (!success) {
		fout << "Error creating directory" << endl;
	}

	ofstream ofile;
	string ofilename = path + "/test_dirichlet.data";
	ofile.open(ofilename);

	for (int i = 0; i < N; ++i) {
		Suffies_MultivariateNormal *suffies = (Suffies_MultivariateNormal*)hyper(generator);
		ofile << suffies->mu << endl;
	}
	ofile.close();

}
