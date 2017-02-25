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

	int D = 2;

	// number of samples to generate
	int N = 100;

	default_random_engine generator;
	
	fout << "Dirichlet distribution" << endl;
	Suffies_Dirichlet suffies_dirichlet;
	suffies_dirichlet.alpha = 1;

	fout << "Normal Inverse Wishart distribution" << endl;
	Suffies_NormalInvWishart *suffies_niw;
	string ofilename = "";
	switch (D) {
		case 1: 
		{
			suffies_niw = new Suffies_NormalInvWishart(1);
			suffies_niw->mu << 6;
			suffies_niw->kappa = 1.0;
			suffies_niw->nu = 4;
			suffies_niw->Lambda << 0.01;

			ofilename = "test_dirichlet.data";
			break;
		}
		case 2: 
		{
			suffies_niw = new Suffies_NormalInvWishart(2);
			suffies_niw->mu << 6, 6;
			fout << "Hyper mu is: " << suffies_niw->mu.transpose() << endl;
			suffies_niw->kappa = 0.001;
			suffies_niw->nu = 4;
			suffies_niw->Lambda << 0.01, 0, 0, 0.01;

			ofilename = "test_dirichlet_2d.data";
			break;
		}
		break;
	}
	normal_inverse_wishart_distribution prior(*suffies_niw);

	dirichlet_process hyper(suffies_dirichlet, prior);

	bool success;

	const string & path = "output/test";
	if (!fs::exists(path)) {
		success = fs::create_directories(path);
		if (!success) {
			fout << "Error creating directory" << endl;
		}
	}

	ofstream ofile;
	ofilename = path + '/' + ofilename;
	ofile.open(ofilename);

	std::vector<Suffies_MultivariateNormal *> history;
	history.clear();
	for (int i = 0; i < N; ++i) {
		Suffies_MultivariateNormal *suffies_mvn = hyper.sample(generator, history);
		history.push_back(suffies_mvn);
	
		multivariate_normal_distribution likelihood(*suffies_mvn);
		Suffies_Unity_MultivariateNormal *suffies = likelihood(generator);

		//ofile << suffies->mu.transpose() << endl;
		ofile << suffies->mu.transpose() << endl;
	}
	ofile.close();

	fout << "Written results to " << ofilename << endl;

	delete suffies_niw;
}
