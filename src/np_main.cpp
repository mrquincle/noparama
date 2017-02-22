#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>

#include <np_mcmc.h>
#include <np_data.h>
#include <np_results.h>

#include <membertrix.h>

#include <statistics/multivariatenormal.h>
#include <statistics/dirichlet.h>
#include <statistics/normalinvwishart.h>

#include <np_neal_algorithm8.h>

#include <pretty_print.hpp>

using namespace std;

int main(int argc, char *argv[]) {
	char _verbosity = 2;
	
	fout << "Welcome to noparama" << endl;

	// Configuration parameters 
	int T = 1000;
	double alpha = 1;
	
	default_random_engine generator(random_device{}()); 

	// Load data
	string datafilename;
	if (argc > 1) {
		datafilename = string(argv[1]);
	}
	else {
		fout << "The datafile is required as argument." << endl;
		fout << "For example:" << endl;
		datafilename = "$HOME/workspace/thesis/dpm/inference/data/many-modal/pattern100_sigma0_1.plain.txt";
		fout << argv[0] << ' ' << datafilename << endl;
		fout << "Or:" << endl;
		datafilename = "datasets/twogaussians.data";
		fout << argv[0] << ' ' << datafilename << endl;
		exit(1);
	}
	fout << "Load file: " << datafilename << endl;

	dataset_t dataset;
	// The data file should have "a b c" on lines, separated by spaces (without quotes, each value of the type double).
	fout << "Read dataset" << endl;
	std::ifstream datafilehandle(datafilename);
	double a, b,c;
	ground_truth_t ground_truth;
	ground_truth.clear();
	int n = 0;
	while (datafilehandle >> a >> b >> c) {
		data_t *data = new data_t(2);
		*data = { a, b };
		dataset.push_back(data);
		ground_truth[n] = (int)c;
		n++;
	}

	int I = 20;
	fout << "Display first " << I << " items of the dataset" << endl;
	for (int i = 0; i < I; ++i) {
		if (i == (int)ground_truth.size()) break;
		fout << "Data: " << *dataset[i] << " with ground truth " << ground_truth[i] << endl;
	}

	// The likelihood function is a multivariate normal distribution (parameters need not be set)
	fout << "Multivariate normal suffies" << endl;
	Suffies_MultivariateNormal suffies_mvn(2);
	suffies_mvn.mu << 0, 0;
	suffies_mvn.sigma << 1, 0, 0, 1;
	fout << "Multivariate normal distribution" << endl;
	multivariate_normal_distribution likelihood(suffies_mvn);

	// The hierarchical prior has an alpha of 1 and the base distribution is handed separately through the prior
	fout << "Dirichlet distribution" << endl;
	Suffies_Dirichlet suffies_dirichlet;
	suffies_dirichlet.alpha = alpha;

	fout << "Normal Inverse Wishart distribution" << endl;
	Suffies_NormalInvWishart suffies_niw(2);
	suffies_niw.mu << 6, 6;
	suffies_niw.kappa = 1.0/500;
	suffies_niw.nu = 4;
	suffies_niw.Lambda << 0.01, 0, 0, 0.01;

	normal_inverse_wishart_distribution prior(suffies_niw);

	dirichlet_distribution hyper(suffies_dirichlet, prior);
	
	fout << "Set up InitClusters object" << endl;
	InitClusters init_clusters(generator, hyper);

	// To update cluster parameters we need prior (hyper parameters) and likelihood, the membership matrix is left 
	// invariant
	fout << "Set up UpdateClusters object" << endl;
	UpdateClusters update_clusters(generator, likelihood, hyper);

	// To update the cluster population we need	to sample new clusters using hyper parameters and adjust existing ones
	// using prior and likelihood
	fout << "Set up UpdateClusterPopulation object" << endl;
	NealAlgorithm8 update_cluster_population(generator, likelihood, hyper);

	// create MCMC object
	fout << "Set up MCMC" << endl;
	MCMC & mcmc = *new MCMC(generator, init_clusters, update_clusters, update_cluster_population);

	fout << "Run MCMC for " << T << " steps" << endl;
	mcmc.run(dataset, T);

	update_cluster_population.printStatistics();

	const membertrix trix = mcmc.getMembershipMatrix();

	// analyse and write out results
	std::chrono::time_point<std::chrono::system_clock> clock;
	clock = std::chrono::system_clock::now();
	std::time_t time = std::chrono::system_clock::to_time_t(clock);
	std::tm tm = *std::localtime(&time);

	std::string workspace = "output/";
	std::stringstream tss; 
	tss << std::put_time(&tm, "%Y%m%d_%H:%M");
	std::string dirname = tss.str(); 
	std::string basename = "results";

	Results results(trix, ground_truth);
	results.write(workspace, dirname, basename);

}


