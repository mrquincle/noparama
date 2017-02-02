#include <iostream>
#include <fstream>

#include <np_mcmc.h>
#include <np_data.h>

#include <statistics/multivariatenormal.h>
#include <statistics/dirichlet.h>
#include <statistics/normalinvwishart.h>

int main() {
	std::cout << "Welcome to noparama" << std::endl;
	std::string datafilename;

	// Load data
	datafilename = "/home/anne/workspace/thesis/dpm/inference/data/many-modal/pattern100_sigma0_1.plain.txt";

	dataset_t dataset;

	// The data file should have "a b" on lines, separated by spaces (without quotes, each value of the type double).
	std::ifstream datafilehandle(datafilename);
	double a, b;
	while (datafilehandle >> a >> b) {
		data_t &data = *new data_t(2);
		data = { a, b };
		dataset.push_back(&data);
	}

	// The likelihood function is a multivariate normal distribution (parameters need not be set)
	Suffies_MultivariateNormal suffies_mvn;
	multivariate_normal_distribution likelihood(suffies_mvn);

	// The hierarchical prior has an alpha of 1 and the base distribution is handed separately through the prior
	Suffies_Dirichlet suffies_dirichlet;
	suffies_dirichlet.alpha = 1;

	// TODO: check how we can add the niw to hyper distribution, probably not via suffies but via distribution itself

	dirichlet_distribution hyper(suffies_dirichlet);

	// To update cluster parameters we need prior (hyper parameters) and likelihood, the membership matrix is left 
	// invariant
	UpdateClusters update_clusters( (distribution_t&)likelihood, (distribution_t&)hyper);

	Suffies_NormalInvWishart suffies_niw;
	suffies_niw.mu << 6, 6;
	suffies_niw.kappa = 1/500;
	suffies_niw.nu = 4;
	suffies_niw.Lambda << 0.01, 0, 0, 0.01;

	normal_inverse_wishart_distribution prior(suffies_niw);

	// To update the cluster population we need	to sample new clusters using hyper parameters and adjust existing ones
	// using prior and likelihood
	UpdateClusterPopulation update_cluster_population(prior, likelihood, hyper);

	// create MCMC object
	MCMC & mcmc = *new MCMC(update_clusters, update_cluster_population);

	mcmc.run(dataset);

}


