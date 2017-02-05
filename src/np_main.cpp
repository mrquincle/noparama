#include <iostream>
#include <fstream>

#include <np_mcmc.h>
#include <np_data.h>

#include <statistics/multivariatenormal.h>
#include <statistics/dirichlet.h>
#include <statistics/normalinvwishart.h>

#include <np_neal_algorithm8.h>

using namespace std;

int main() {
	cout << "Welcome to noparama" << endl;
	string datafilename;

	// Load data
	datafilename = "/home/anne/workspace/thesis/dpm/inference/data/many-modal/pattern100_sigma0_1.plain.txt";

	dataset_t dataset;

	// The data file should have "a b" on lines, separated by spaces (without quotes, each value of the type double).
	cout << "Read dataset" << endl;
	std::ifstream datafilehandle(datafilename);
	double a, b;
	while (datafilehandle >> a >> b) {
		data_t &data = *new data_t(2);
		data = { a, b };
		dataset.push_back(&data);
	}

	// The likelihood function is a multivariate normal distribution (parameters need not be set)
	cout << "Multivariate normal suffies" << endl;
	Suffies_MultivariateNormal suffies_mvn(2);
	suffies_mvn.mu << 0, 0;
	suffies_mvn.sigma << 1, 0, 0, 1;
	cout << "Multivariate normal distribution" << endl;
	multivariate_normal_distribution likelihood(suffies_mvn);

	// The hierarchical prior has an alpha of 1 and the base distribution is handed separately through the prior
	cout << "Dirichlet distribution" << endl;
	Suffies_Dirichlet suffies_dirichlet;
	suffies_dirichlet.alpha = 1;

	cout << "Normal Inverse Wishart distribution" << endl;
	Suffies_NormalInvWishart suffies_niw(2);
	suffies_niw.mu << 6, 6;
	suffies_niw.kappa = 1/500;
	suffies_niw.nu = 4;
	suffies_niw.Lambda << 0.01, 0, 0, 0.01;

	normal_inverse_wishart_distribution prior(suffies_niw);

	dirichlet_distribution hyper(suffies_dirichlet, prior);

	// To update cluster parameters we need prior (hyper parameters) and likelihood, the membership matrix is left 
	// invariant
	cout << "Set up UpdateClusters object" << endl;
	UpdateClusters update_clusters( (distribution_t&)likelihood, (distribution_t&)hyper);

	// To update the cluster population we need	to sample new clusters using hyper parameters and adjust existing ones
	// using prior and likelihood
	cout << "Set up UpdateClusterPopulation object" << endl;
	NealAlgorithm8 update_cluster_population(likelihood, hyper);

	// create MCMC object
	cout << "Set up MCMC" << endl;
	MCMC & mcmc = *new MCMC(update_clusters, update_cluster_population);

	cout << "Run MCMC" << endl;
	mcmc.run(dataset);

}


