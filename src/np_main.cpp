#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <unordered_map>

// Include getopt
#include <unistd.h>
#include <experimental/filesystem>

#include <np_mcmc.h>
#include <np_data.h>
#include <np_results.h>

#include <membertrix.h>

#include <statistics/scalarnoise_multivariatenormal.h>
#include <statistics/multivariatenormal.h>
#include <statistics/dirichlet.h>
#include <statistics/normalinvwishart.h>
#include <statistics/normalinvgamma.h>

#include <data-driven/set_compare.h>

//#include <np_neal_algorithm2.h>
#include <np_neal_algorithm8.h>
#include <np_jain_neal_algorithm.h>
#include <np_triadic_algorithm.h>

#include <pretty_print.hpp>

namespace fs = std::experimental::filesystem;

enum algorithm_t { /*algorithm2, */ algorithm8, jain_neal_split, triadic };

void disp_help(std::string appname) {
	std::string datafilename;
	std::cout << "noparama [version 0.1.72]" << std::endl;
	std::cout << std::endl;
	std::cout << "Usage: " << std::endl;
	std::cout << "  noparama -d datafile -a algorithm8|jain_neal_split|triadic  [-T ticks] [-c regression|clustering|angular]" << std::endl;
	std::cout << "For example:" << std::endl;
	datafilename = "datasets/twogaussians.data";
	std::cout << "  " << appname << " -d " << datafilename << " -a triadic -T 5000 -c regression" << std::endl;
}

/**
 * With the dataset twogaussians it is important to realize that this did not originate from a Dirichlet process. Even
 * if alpha is chosen to be very low (say 0.0001), there will be multiple clusters generated, not just two. The purity
 * will be quite high (very few misclassifications) but the specificity is low (identification of multiple clusters
 * where there is only one).
 */
int main(int argc, char *argv[]) {
	char _verbosity = Debug;
	
	fout << "Welcome to noparama" << std::endl;

	// ---------------------------------------------------------------------------------------------------------------
	// Configuration parameters 
	// ---------------------------------------------------------------------------------------------------------------
	double alpha = 1;
	// ---------------------------------------------------------------------------------------------------------------

	fout << "The Dirichlet Process is run with alpha=" << alpha << std::endl;

	// limit number of data items
	bool limit = false;
	int N = 10;
	if (limit) {
		fout << "Using limited dataset with only " << N << " items " << std::endl;
	} else {
		fout << "Using full dataset" << std::endl;
	}

	std::default_random_engine generator(std::random_device{}()); 

	// Configuration std::strings that can be filled from CLI arguments
	std::string datafilename;
	int token;
	std::unordered_map<int, std::string> flags;
	while( (token = getopt(argc, argv, "d:c:a:T:h?")) != EOF)
	{
		switch (token)
		{
			case 'd':
				flags['d'] = std::string(optarg);
				break;
			case 'c':
				flags['c'] = std::string(optarg);
				break;
			case 'a':
				flags['a'] = std::string(optarg);
				break;
			case 'T':
				flags['T'] = std::string(optarg);
				break;
			case 'h': case '?':
				disp_help(std::string(argv[0]));
				exit(1);
				break;
		}
	}

	// required, so show help and exit
	std::string algorithm_str = flags['a'];
	algorithm_t algorithm;
	if (flags.find('d') == flags.end() || flags.find('a') == flags.end()) {
		disp_help(std::string(argv[0]));
		exit(1);
	} else {
		datafilename = flags['d'];

		/*if (algorithm_str == "algorithm2") {
			algorithm = algorithm2;
		} else */
		if (algorithm_str == "algorithm8") {
			algorithm = algorithm8;
		} else if (algorithm_str == "jain_neal_split") {
			algorithm = jain_neal_split;
		} else if (algorithm_str == "triadic") {
			algorithm = triadic;
		} else {
			std::cerr << "Unknown algorithm: " << algorithm_str << std::endl;
			exit(1);
		}
	}
	
	// optional T flag
	int T = 2000;
	if (flags.find('T') != flags.end()) {
		T = std::stoi(flags['T']);
	}

	// optional clustering/regression flag
	representation_mode_t representation_mode = regression_mode; 
	if (flags.find('c') != flags.end()) {
		std::string configuration = flags['c'];
		if (configuration == "regression") {
			fout << "Regression mode" << std::endl;
			representation_mode = regression_mode;
		} else if (configuration == "clustering") {
			fout << "Clustering mode" << std::endl;
			representation_mode = clustering_mode;
		} else if (configuration == "angular") {
			fout << "Angular mode (form of regression)" << std::endl;
			representation_mode = angular_mode;
		} else {
			fout << "Picked regression mode by default." << std::endl;
		}
	}

	// check if output directory exists and bail out if this is the case
	std::string justdatafilename = fs::path(datafilename).filename();
	std::string workspace = "output/" + algorithm_str + '/' + justdatafilename + '/';
	fout << "Have workspace: " << workspace << std::endl;
	if (fs::exists(workspace)) {
		std::cerr << "Directory already exists. Drop out. We don't want to overwrite it or do double work." << std::endl;
		exit(106);
	}
	
	fout << "Run MCMC sampler for " << T << " steps " << std::endl;

	fout << "Load file: " << datafilename << std::endl;

	dataset_t dataset;
	// The data file should have "a b c" on lines, separated by spaces (without quotes, each value of the type double).
	fout << "Read dataset" << std::endl;
	std::ifstream datafilehandle(datafilename);
	double a, b, c;
	ground_truth_t ground_truth;
	ground_truth.clear();
	int n = 0;
	while (datafilehandle >> a >> b >> c) {
		switch(representation_mode) {
		case regression_mode:
			{
				data_t *data = new data_t(3);
				// prepend vector with constant for regression
				*data = { 1, a, b };
				dataset.push_back(data);
				break;
			}
		case angular_mode: case clustering_mode: 
			{
				data_t *data = new data_t(2);
				*data = { a, b };
				dataset.push_back(data);
				break;
			}
		}
		ground_truth[n] = (int)c;
		n++;
		if (limit && (n == N)) break;
	}

	if (n == 0) {
		fout << "No data found... Check the file or the contents of the file." << std::endl;
		exit(7);
	}

	int I = 5;
	fout << "Display first " << I << " items of the dataset" << std::endl;
	for (int i = 0; i < I; ++i) {
		if (i == (int)ground_truth.size()) break;
		fout << "Data: " << *dataset[i] << " with ground truth " << ground_truth[i] << std::endl;
	}

	// The likelihood function only is required w.r.t. type, parameters are normally set later
	distribution_t *likelihood;
	if (representation_mode == regression_mode || representation_mode == angular_mode) {
		fout << "Multivariate normal suffies (sigma scalar)" << std::endl;
		// note that we need to allocate the suffies or else it won't exists after the call to new multivar..
		Suffies_ScalarNoise_MultivariateNormal * suffies_mvn = new Suffies_ScalarNoise_MultivariateNormal(2);
		suffies_mvn->mu << 0, 0;
		suffies_mvn->sigma = 1;
		likelihood = new scalarnoise_multivariate_normal_distribution(*suffies_mvn, representation_mode);
	} else {
		fout << "Multivariate normal suffies" << std::endl;
		Suffies_MultivariateNormal * suffies_mvn = new Suffies_MultivariateNormal(2);
		suffies_mvn->mu << 0, 0;
		suffies_mvn->sigma << 1, 0, 0, 1; // get rid of this
		likelihood = new multivariate_normal_distribution(*suffies_mvn);
	}
		
	fout << "likelihood still: " << likelihood->getSuffies() << std::endl;

	// The hierarchical prior has an alpha of 1 and the base distribution is handed separately through the prior
	fout << "Dirichlet distribution" << std::endl;
	Suffies_Dirichlet suffies_dirichlet;
	suffies_dirichlet.alpha = alpha;
	distribution_t *prior;

	if (representation_mode == regression_mode || representation_mode == angular_mode) {
		fout << "Normal Inverse Gamma distribution" << std::endl;
		Suffies_NormalInvGamma * suffies_nig = new Suffies_NormalInvGamma(2);
		suffies_nig->mu << 0, 0;
		suffies_nig->alpha = 10;
		suffies_nig->beta = 0.1;
		suffies_nig->Lambda << 0.01, 0, 0, 0.01;
		prior = new normal_inverse_gamma_distribution(*suffies_nig);
	} else {
		fout << "Normal Inverse Wishart distribution" << std::endl;
		Suffies_NormalInvWishart * suffies_niw = new Suffies_NormalInvWishart(2);
		suffies_niw->mu << 6, 6;
		suffies_niw->kappa = 1.0/500;
		suffies_niw->nu = 4;
		suffies_niw->Lambda << 0.01, 0, 0, 0.01;
		prior = new normal_inverse_wishart_distribution(*suffies_niw);
	}
	fout << "Create Dirichlet Process" << std::endl;
	dirichlet_process hyper(suffies_dirichlet, *prior);
	
	fout << "Set up InitClusters object" << std::endl;
	InitClusters init_clusters(generator, hyper);

	// To update cluster parameters we need prior (hyper parameters) and likelihood, the membership matrix is left 
	// invariant
	fout << "Set up UpdateClusters object" << std::endl;
	UpdateClusters update_clusters(generator, *likelihood, hyper);

	// To update the cluster population we need	to sample new clusters using hyper parameters and adjust existing ones
	// using prior and likelihood
	fout << "Set up UpdateClusterPopulation object" << std::endl;

	int subset_count;
	UpdateClusterPopulation *update_cluster_population;
	switch(algorithm) {
/*		case algorithm2: 
			{
				fout << "We will be using algorithm 2 by Neal" << std::endl;
				update_cluster_population = new NealAlgorithm2(generator, *likelihood, hyper);
				subset_count = 1;
				break;
			}
*/
		case algorithm8: 
			{
				fout << "We will be using algorithm 8 by Neal" << std::endl;
				update_cluster_population = new NealAlgorithm8(generator, *likelihood, hyper);
				subset_count = 1;
				break;
			}
		case jain_neal_split:
			{
				fout << "We will be using the Jain-Neal algorithm" << std::endl;
				update_cluster_population = new JainNealAlgorithm(generator, *likelihood, hyper);
				subset_count = 2;
				break;
			}
		case triadic:
			{
				fout << "We will be using the Triadic algorithm" << std::endl;
				update_cluster_population = new TriadicAlgorithm(generator, *likelihood, hyper);
//#define TEST_WITH_TWO
#ifdef TEST_WITH_TWO
				subset_count = 2;
#else
				subset_count = 3;
#endif
				break;
			}
	}

	// create MCMC object
	fout << "Set up MCMC" << std::endl;
	MCMC & mcmc = *new MCMC(generator, init_clusters, update_clusters, *update_cluster_population, subset_count, *likelihood);

	fout << "Run MCMC for " << T << " steps" << std::endl;
	mcmc.run(dataset, T);

	fout << "Print statistics" << std::endl;
	update_cluster_population->printStatistics();
	
	// analyse and write out results
	std::chrono::time_point<std::chrono::system_clock> clock;
	clock = std::chrono::system_clock::now();
	std::time_t time = std::chrono::system_clock::to_time_t(clock);
	std::tm tm = *std::localtime(&time);

	std::stringstream tss; 
	tss << std::put_time(&tm, "%Y%m%d_%H:%M");
	std::string dirname = tss.str(); 
	std::string basename = "snapshot";
	
	fout << "Last item" << std::endl;
	const membertrix trix_snapshot = mcmc.getMembershipMatrix();
	Results results_snapshot(trix_snapshot, ground_truth);
	results_snapshot.write(workspace, dirname, basename);
	
	basename = "results";
	
	fout << "Write results" << std::endl;
	const membertrix trix = mcmc.getMaxLikelihoodMatrix();
	Results results(trix, ground_truth);
	results.write(workspace, dirname, basename);
	
	fout << "Memory deallocation" << std::endl;
	delete &mcmc;
	delete likelihood;
	delete prior;
	for (int i = 0; i < N; ++i) {
		delete dataset[i];
	}
	delete update_cluster_population;
}


