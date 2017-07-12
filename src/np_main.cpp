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

#include <np_neal_algorithm8.h>
#include <np_jain_neal_algorithm.h>
#include <np_triadic_algorithm.h>

#include <pretty_print.hpp>

using namespace std;

enum algorithm_t { algorithm8, jain_neal_split, triadic };

void disp_help(std::string appname) {
	std::string datafilename;
	cout << "noparama [version 0.1.72]" << endl;
	cout << endl;
	cout << "Usage: " << endl;
	cout << "  noparama -d datafile -a algorithm8|jain_neal_split|triadic  [-T ticks] [-c regression|clustering]" << endl;
	cout << "For example:" << endl;
	datafilename = "datasets/twogaussians.data";
	cout << "  " << appname << " -d " << datafilename << " -a triadic -T 5000 -c regression" << endl;
}

/**
 * With the dataset twogaussians it is important to realize that this did not originate from a Dirichlet process. Even
 * if alpha is chosen to be very low (say 0.0001), there will be multiple clusters generated, not just two. The purity
 * will be quite high (very few misclassifications) but the specificity is low (identification of multiple clusters
 * where there is only one).
 */
int main(int argc, char *argv[]) {
	char _verbosity = Debug;
	
	fout << "Welcome to noparama" << endl;

	// ---------------------------------------------------------------------------------------------------------------
	// Configuration parameters 
	// ---------------------------------------------------------------------------------------------------------------
	double alpha = 1;
	// ---------------------------------------------------------------------------------------------------------------


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
	int token;
	unordered_map<int, string> flags;
	while( (token = getopt(argc, argv, "d:c:a:T:h?")) != EOF)
	{
		switch (token)
		{
			case 'd':
				flags['d'] = string(optarg);
				break;
			case 'c':
				flags['c'] = string(optarg);
				break;
			case 'a':
				flags['a'] = string(optarg);
				break;
			case 'T':
				flags['T'] = string(optarg);
				break;
			case 'h': case '?':
				disp_help(std::string(argv[0]));
				exit(1);
				break;
		}
	}

	// required, so show help and exit
	algorithm_t algorithm;
	if (flags.find('d') == flags.end() || flags.find('a') == flags.end()) {
		disp_help(std::string(argv[0]));
		exit(1);
	} else {
		datafilename = flags['d'];

		std::string algorithm_str = flags['a'];
		if (algorithm_str == "algorithm8") {
			algorithm = algorithm8;
		} else if (algorithm_str == "jain_neal_split") {
			algorithm = jain_neal_split;
		} else if (algorithm_str == "triadic") {
			algorithm = triadic;
		} else {
			cerr << "Unknown algorithm: " << algorithm_str << endl;
			exit(1);
		}
	}
	
	// optional T flag
	int T = 2000;
	if (flags.find('T') != flags.end()) {
		T = std::stoi(flags['T']);
	}

	// optional clustering/regression flag
	bool regression = false;
	if (flags.find('c') != flags.end()) {
		std::string configuration = flags['c'];
		if (configuration == "regression") {
			regression = true;
		}
	}
	
	fout << "Run MCMC sampler for " << T << " steps " << endl;

	if (regression) {
		fout << "Regression mode" << endl;
	} else {
		fout << "Clustering mode" << endl;
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
	UpdateClusterPopulation *update_cluster_population;
	switch(algorithm) {
		case algorithm8: 
			{
				fout << "We will be using algorithm 8 by Neal" << endl;
				update_cluster_population = new NealAlgorithm8(generator, *likelihood, hyper);
				subset_count = 1;
				break;
			}
		case jain_neal_split:
			{
				fout << "We will be using the Jain-Neal algorithm" << endl;
				update_cluster_population = new JainNealAlgorithm(generator, *likelihood, hyper);
				subset_count = 2;
				break;
			}
		case triadic:
			{
				fout << "We will be using the Triadic algorithm" << endl;
				update_cluster_population = new TriadicAlgorithm(generator, *likelihood, hyper);
#define TEST_WITH_TWO
#ifdef TEST_WITH_TWO
				subset_count = 2;
#else
				subset_count = 3;
#endif
				break;
			}
	}

	// create MCMC object
	fout << "Set up MCMC" << endl;
	MCMC & mcmc = *new MCMC(generator, init_clusters, update_clusters, *update_cluster_population, subset_count, *likelihood);

	fout << "Run MCMC for " << T << " steps" << endl;
	mcmc.run(dataset, T);

	fout << "Print statistics" << endl;
	update_cluster_population->printStatistics();

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
	delete update_cluster_population;
}


