#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <unordered_map>

// Include getopt
#include <unistd.h>

#include <np_mcmc.h>
#include <np_data.h>
#include <np_results.h>

#include <membertrix.h>

#include <statistics/scalarnoise_multivariatenormal.h>
#include <statistics/multivariatenormal.h>
#include <statistics/dirichlet.h>
#include <statistics/normalinvwishart.h>
#include <statistics/normalinvgamma.h>

#define ALGORITHM 9

#if ALGORITHM==8
#include <np_neal_algorithm8.h>
#else
#include <np_jain_neal_algorithm.h>
#endif

#include <pretty_print.hpp>

using namespace std;

/**
 * With the dataset twogaussians it is important to realize that this did not originate from a Dirichlet process. Even
 * if alpha is chosen to be very low (say 0.0001), there will be multiple clusters generated, not just two. The purity
 * will be quite high (very few misclassifications) but the specificity is low (identification of multiple clusters
 * where there is only one).
 */
int main(int argc, char *argv[]) {
	char _verbosity = Debug;
	
	fout << "Welcome to noparama" << endl;

	// Configuration parameters 
	int T = 1000;
	fout << "Run MCMC sampler for " << T << " steps " << endl;

	double alpha = 1;
	fout << "The Dirichlet Process is run with alpha=" << alpha << endl;

	// limit number of data items
	bool limit = false;
	int N = 10;
	if (limit) {
		fout << "Using limited dataset with only " << N << " items " << endl;
	} else {
		fout << "Using full dataset" << endl;
	}

	default_random_engine generator(random_device{}()); 

	// Configuration strings that can be filled from CLI arguments
	string datafilename;
	string configuration;

	int token;
	unordered_map<int, string> flags;
	while( (token = getopt(argc, argv, "d:c:")) != EOF)
	{
		switch (token)
		{
			case 'd':
				flags['d'] = string(optarg);
				break;
			case 'c':
				flags['c'] = string(optarg);
				break;
		}
	}

	if (flags.find('d') == flags.end()) {
		// required, so show help and exit
		cout << "noparama [version 0.0.1]" << endl;
		cout << endl << endl;
		cout << "Usage: " << endl;
		fout << "  noparama [arguments]" << endl;
		datafilename = "$HOME/workspace/thesis/dpm/inference/data/many-modal/pattern100_sigma0_1.plain.txt";
		fout << argv[0] << " -d " << datafilename << endl;
		fout << "Or:" << endl;
		datafilename = "datasets/twogaussians.data";
		fout << argv[0] << " -d " << datafilename << endl;
		exit(1);
	} else {
		datafilename = flags['d'];
	}

	bool regression = false;
	if (flags.find('c') == flags.end()) {
		// no configuration parameter, but optional anyway
			fout << "Clustering mode" << endl;
	} else {
		configuration = flags['c'];
		if (configuration == "regression") {
			regression = true;
			fout << "Regression mode" << endl;
		}
	}

	fout << "Load file: " << datafilename << endl;

	dataset_t dataset;
	// The data file should have "a b c" on lines, separated by spaces (without quotes, each value of the type double).
	fout << "Read dataset" << endl;
	std::ifstream datafilehandle(datafilename);
	double a, b, c;
	ground_truth_t ground_truth;
	ground_truth.clear();
	int n = 0;
	while (datafilehandle >> a >> b >> c) {
		if (regression) {
			// prepend vector with constant for regression
			data_t *data = new data_t(3);
			*data = { 1, a, b };
			dataset.push_back(data);
		} else {
			data_t *data = new data_t(2);
			*data = { a, b };
			dataset.push_back(data);
		}
		ground_truth[n] = (int)c;
		n++;
		if (limit && (n == N)) break;
	}

	if (n == 0) {
		fout << "No data found... Check the file or the contents of the file." << endl;
		exit(7);
	}

	int I = 5;
	fout << "Display first " << I << " items of the dataset" << endl;
	for (int i = 0; i < I; ++i) {
		if (i == (int)ground_truth.size()) break;
		fout << "Data: " << *dataset[i] << " with ground truth " << ground_truth[i] << endl;
	}

	// The likelihood function only is required w.r.t. type, parameters are normally set later
	distribution_t *likelihood;
	if (regression) {
		fout << "Multivariate normal suffies (sigma scalar)" << endl;
		// note that we need to allocate the suffies or else it won't exists after the call to new multivar..
		//Suffies_ScalarNoise_MultivariateNormal * suffies_mvn = new Suffies_ScalarNoise_MultivariateNormal(2);
		Suffies_ScalarNoise_MultivariateNormal * suffies_mvn = new Suffies_ScalarNoise_MultivariateNormal(2);
		suffies_mvn->mu << 0, 0;
		suffies_mvn->sigma = 1;
		likelihood = new scalarnoise_multivariate_normal_distribution(*suffies_mvn, true);
	} else {
		fout << "Multivariate normal suffies" << endl;
		Suffies_MultivariateNormal * suffies_mvn = new Suffies_MultivariateNormal(2);
		suffies_mvn->mu << 0, 0;
		suffies_mvn->sigma << 1, 0, 0, 1; // get rid of this
		likelihood = new multivariate_normal_distribution(*suffies_mvn);
	}
		
	fout << "likelihood still: " << likelihood->getSuffies() << endl;

	// The hierarchical prior has an alpha of 1 and the base distribution is handed separately through the prior
	fout << "Dirichlet distribution" << endl;
	Suffies_Dirichlet suffies_dirichlet;
	suffies_dirichlet.alpha = alpha;
	distribution_t *prior;

	if (regression) {
		fout << "Normal Inverse Gamma distribution" << endl;
		Suffies_NormalInvGamma * suffies_nig = new Suffies_NormalInvGamma(2);
		suffies_nig->mu << 0, 0;
		suffies_nig->alpha = 10;
		suffies_nig->beta = 0.1;
		suffies_nig->Lambda << 0.01, 0, 0, 0.01;
		prior = new normal_inverse_gamma_distribution(*suffies_nig);
	} else {
		fout << "Normal Inverse Wishart distribution" << endl;
		Suffies_NormalInvWishart * suffies_niw = new Suffies_NormalInvWishart(2);
		suffies_niw->mu << 6, 6;
		suffies_niw->kappa = 1.0/500;
		suffies_niw->nu = 4;
		suffies_niw->Lambda << 0.01, 0, 0, 0.01;
		prior = new normal_inverse_wishart_distribution(*suffies_niw);
	}
	fout << "Create Dirichlet Process" << endl;
	dirichlet_process hyper(suffies_dirichlet, *prior);
	
	fout << "Set up InitClusters object" << endl;
	InitClusters init_clusters(generator, hyper);

	// To update cluster parameters we need prior (hyper parameters) and likelihood, the membership matrix is left 
	// invariant
	fout << "Set up UpdateClusters object" << endl;
	UpdateClusters update_clusters(generator, *likelihood, hyper);

	// To update the cluster population we need	to sample new clusters using hyper parameters and adjust existing ones
	// using prior and likelihood
	fout << "Set up UpdateClusterPopulation object" << endl;

	int subset_count;
#if ALGORITHM==8
	fout << "We will be using algorithm 8 by Neal" << endl;
	NealAlgorithm8 update_cluster_population(generator, *likelihood, hyper);
	subset_count = 1;
#else
	fout << "We will be using the Jain-Neal algorithm" << endl;
	JainNealAlgorithm update_cluster_population(generator, *likelihood, hyper);
	subset_count = 2;
#endif
	// create MCMC object
	fout << "Set up MCMC" << endl;
	MCMC & mcmc = *new MCMC(generator, init_clusters, update_clusters, update_cluster_population, subset_count);

	fout << "Run MCMC for " << T << " steps" << endl;
	mcmc.run(dataset, T);

	fout << "Print statistics" << endl;
	update_cluster_population.printStatistics();

	fout << "Write results" << endl;
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
	
	fout << "Memory deallocation" << endl;

	delete &mcmc;
	delete likelihood;
	delete prior;
	for (int i = 0; i < N; ++i) {
		delete dataset[i];
	}
}


